from flask import Flask, request, jsonify, render_template
from Shor_Code import ShorCode

app = Flask(__name__)

@app.route('/')
def index():
    # serves the visualiser page
    return render_template('index.html')

@app.route('/run_circuit', methods=['POST'])
def run_circuit():
    data = request.json
    x_error = data.get('bitError')    # e.g., {'block': 1, 'idx': 1}
    z_error = data.get('phaseError')  # e.g., 2
    
    # run the Shor code with the chosen errors
    code = ShorCode()
    results = code.Execute(x_pos=x_error, z_pos=z_error)
    
    # send the results back to the page
    return jsonify({
        "status": "success",
        "simulation_results": results
    })

if __name__ == '__main__':
    app.run(debug=True, port=5000)
