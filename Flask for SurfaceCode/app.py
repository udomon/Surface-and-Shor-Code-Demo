"""
Flask backend for the Surface Code demo.

All the decoding (layout, stabilisers, error injection, syndromes, matching
and corrections) is done by the classes in Surface_Code.py, the JavaScript
frontend only draws what this backend sends back.
"""

import uuid
from flask import Flask, request, jsonify, render_template
import networkx as nx

from Surface_Code import Error_Correction, Surface_Code_Graph

app = Flask(__name__)

# stores the Error_Correction object, plan and current step for each session
SESSIONS = {}


def coord_list(d):
    """[i, j] -> (i, j)"""
    return tuple(d)


def serialize_coords(qubit_dict):
    return [[i, j] for (i, j) in qubit_dict.keys()]


def serialize_stabilisers(stab_dict):
    return [{"anc": [i, j], "neighbours": [[ni, nj] for (ni, nj) in nbs]}
            for (i, j), nbs in stab_dict.items()]


def serialize_syndrome(syn_dict):
    return {f"{i},{j}": v for (i, j), v in syn_dict.items()}


def serialize_error_set(err_set):
    return [[i, j] for (i, j) in err_set]


def layout_payload(ec):
    return {
        "d": ec.d,
        "size": 2 * ec.d - 1,
        "data_qubits": serialize_coords(ec.data_q),
        "x_anc": serialize_coords(ec.x_anc),
        "z_anc": serialize_coords(ec.z_anc),
        "x_stabilisers": serialize_stabilisers(ec.x_stabilisers),
        "z_stabilisers": serialize_stabilisers(ec.z_stabilisers),
    }


def state_payload(ec):
    return {
        "x_errors": serialize_error_set(ec.x_errors),
        "z_errors": serialize_error_set(ec.z_errors),
        "y_errors": serialize_error_set(ec.y_errors),
        "x_syndrome": serialize_syndrome(ec.x_syndrome),
        "z_syndrome": serialize_syndrome(ec.z_syndrome),
    }


def compute_pair_path(ec, value_1, value_2):
    """
    Same path calculation as Error_Correction._flip_bits but for one pair,
    and without changing ec._touched.
    """
    path_cells = set()

    left, right = ec._sort_pairs(value_1, value_2, 'j')
    row = left[0]
    for j in range(left[1], right[1] + 1):
        coord = (row, j)
        if coord in ec.data_q:
            path_cells.add(coord)

    top, bottom = ec._sort_pairs(left, right, 'i')
    col = right[1]
    for i in range(top[0], bottom[0] + 1):
        coord = (i, col)
        if coord in ec.data_q:
            path_cells.add(coord)

    return path_cells


def build_plan_for_type(ec, syndrome_type, syndrome_dict):
    """
    Builds the matching graph the same way as
    Error_Correction.minimum_weight_matching, but returns a list of steps
    instead of applying all the corrections at once. Doesn't change the
    errors in ec, that happens one step at a time later.
    """
    errors = [syn for syn, val in syndrome_dict.items() if val == 1]

    G = nx.Graph()
    size = 2 * ec.d - 1

    boundary_ids = []
    for error_coord in errors:
        G.add_node(error_coord)
        boundary_id = ('boundary', error_coord)
        boundary_ids.append(boundary_id)

        i, j = error_coord
        if syndrome_type == 'X':
            boundary_distance = min(j, size - 1 - j) + 0.5
        else:
            boundary_distance = min(i, size - 1 - i) + 0.5
        G.add_edge(error_coord, boundary_id, weight=-boundary_distance)

    for i in range(len(errors)):
        for j in range(i + 1, len(errors)):
            G.add_edge(errors[i], errors[j], weight=-ec._manhattan_distance(errors[i], errors[j]))

    for i in range(len(boundary_ids)):
        for j in range(i + 1, len(boundary_ids)):
            G.add_edge(boundary_ids[i], boundary_ids[j], weight=0)

    matches = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)

    steps = []
    # sorted so the steps always come out in the same order
    for pair in sorted(matches, key=lambda p: str(p)):
        value_1, value_2 = pair
        i_1, j_1 = value_1
        i_2, j_2 = value_2

        if i_1 == 'boundary' and i_2 == 'boundary':
            continue  # unused boundary-boundary pairing, no physical effect

        raw_from, raw_to, boundary = None, None, False
        resolved_1, resolved_2 = value_1, value_2

        if i_1 == 'boundary':
            raw_from, raw_to = value_2, None
            boundary = True
            resolved_1 = ec._get_boundary_coord(value_2, syndrome_type)
        elif i_2 == 'boundary':
            raw_from, raw_to = value_1, None
            boundary = True
            resolved_2 = ec._get_boundary_coord(value_1, syndrome_type)
        else:
            raw_from, raw_to = value_1, value_2

        path = compute_pair_path(ec, resolved_1, resolved_2)
        if not path:
            continue

        steps.append({
            "type": syndrome_type,
            "raw_from": list(raw_from),
            "raw_to": list(raw_to) if raw_to else None,
            "boundary": boundary,
            "boundary_point": list(ec._get_boundary_coord(raw_from, syndrome_type)) if boundary else None,
            "path": [list(c) for c in sorted(path)],
        })

    return steps


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/api/layout', methods=['POST'])
def api_layout():
    body = request.get_json(force=True)
    d = int(body.get('d', 3))
    if d < 2:
        return jsonify({"error": "d must be >= 2"}), 400
    graph = Surface_Code_Graph(d)
    return jsonify(layout_payload(graph))


@app.route('/api/session/start', methods=['POST'])
def api_session_start():
    body = request.get_json(force=True)
    d = int(body.get('d', 3))
    x_errors = [coord_list(c) for c in body.get('x_errors', [])]
    z_errors = [coord_list(c) for c in body.get('z_errors', [])]
    resync = body.get('resync', 'classical')  # 'classical' for snappy editing clicks, 'quantum' before correcting

    if d < 2:
        return jsonify({"error": "d must be >= 2"}), 400

    ec = Error_Correction(d)

    for coord in x_errors + z_errors:
        if coord not in ec.data_q:
            return jsonify({"error": f"{coord} is not a data qubit, error cannot be inserted"}), 400

    ec.inject_errors(x_pos=x_errors, z_pos=z_errors)
    if resync == 'quantum':
        # recalculate syndromes with the quantum circuit, same as Execute()
        ec.x_syndrome, ec.z_syndrome = ec.calculate_syndromes_from_circuit()

    session_id = str(uuid.uuid4())
    SESSIONS[session_id] = {"ec": ec, "plan": [], "step_index": -1}

    resp = {"session_id": session_id}
    resp.update(layout_payload(ec))
    resp.update(state_payload(ec))
    return jsonify(resp)


@app.route('/api/session/<session_id>/plan', methods=['POST'])
def api_session_plan(session_id):
    session = SESSIONS.get(session_id)
    if not session:
        return jsonify({"error": "unknown session"}), 404
    ec = session["ec"]

    # both plans use the syndromes from before any correction, same as Execute()
    x_plan = build_plan_for_type(ec, 'X', ec.x_syndrome)
    z_plan = build_plan_for_type(ec, 'Z', ec.z_syndrome)

    plan = x_plan + z_plan
    session["plan"] = plan
    session["step_index"] = -1

    return jsonify({
        "x_plan": x_plan,
        "z_plan": z_plan,
        "total_steps": len(plan),
    })


@app.route('/api/session/<session_id>/step', methods=['POST'])
def api_session_step(session_id):
    session = SESSIONS.get(session_id)
    if not session:
        return jsonify({"error": "unknown session"}), 404
    ec = session["ec"]
    plan = session["plan"]

    session["step_index"] += 1
    idx = session["step_index"]
    if idx >= len(plan):
        return jsonify({"error": "no more steps"}), 400

    step = plan[idx]
    for cell in step["path"]:
        coord = tuple(cell)
        ec._apply_correction(coord, step["type"])

    # recalculate syndromes classically so the UI updates between steps
    ec._calculate_syndromes()

    resp = {
        "step_index": idx,
        "total_steps": len(plan),
        "done": (idx + 1) >= len(plan),
        "type": step["type"],
        "raw_from": step["raw_from"],
        "raw_to": step["raw_to"],
        "boundary": step["boundary"],
        "boundary_point": step["boundary_point"],
        "path": step["path"],
    }
    resp.update(state_payload(ec))
    return jsonify(resp)


@app.route('/api/session/<session_id>/finalize', methods=['POST'])
def api_session_finalize(session_id):
    session = SESSIONS.get(session_id)
    if not session:
        return jsonify({"error": "unknown session"}), 404
    ec = session["ec"]

    # final syndrome check with the quantum circuit, same as the end of Execute()
    ec.x_syndrome, ec.z_syndrome = ec.calculate_syndromes_from_circuit()

    return jsonify(state_payload(ec))


@app.route('/api/session/<session_id>/reset', methods=['POST'])
def api_session_reset(session_id):
    SESSIONS.pop(session_id, None)
    return jsonify({"ok": True})


if __name__ == '__main__':
    app.run(debug=True, port=5000)
