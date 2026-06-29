# app.py - simple two-page web shell for asteroid-hunter
import os
import shutil
from flask import Flask, render_template, request, url_for

app = Flask(__name__)

BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
UPLOAD_DIR = os.path.join(BASE_DIR, 'static', 'uploads')
os.makedirs(UPLOAD_DIR, exist_ok=True)

PROJECT_ROOT = os.path.dirname(BASE_DIR)
EXAMPLE_RESULT = os.path.join(PROJECT_ROOT, 'outputs', 'set203_tracks.png')

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/detect', methods=['GET', 'POST'])
def detect():
    result_image = None
    message = None
    if request.method == 'POST':
        uploaded = request.files.get('image')
        if uploaded and uploaded.filename:
            uploaded.save(os.path.join(UPLOAD_DIR, uploaded.filename))
        example_dest = os.path.join(BASE_DIR, 'static', 'example_result.png')
        if os.path.exists(EXAMPLE_RESULT):
            shutil.copy(EXAMPLE_RESULT, example_dest)
            result_image = url_for('static', filename='example_result.png')
            message = ("This is an example detection on the set203 dataset. "
                       "Running detection on your own upload is coming soon.")
        else:
            message = ("Couldn't find the example result image. Run the pipeline "
                       "once (python scripts/exp_set203_pipeline.py) to generate it.")
    return render_template('detect.html', result_image=result_image, message=message)

if __name__ == '__main__':
    app.run(debug=True, port=5000)
