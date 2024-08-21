#!/usr/bin/env python

import os
import subprocess
import argparse
import json
import urllib.request
import hashlib
import sys

def download_models(download_path):
    models_path = os.path.join(sys.prefix, 'bin', 'models.json')
    print(f"Looking for models.json at: {models_path}")  # Debugging line
    if not os.path.exists(models_path):
        raise FileNotFoundError(f"models.json not found at {models_path}")
    
    with open(models_path, 'r') as file:
        models = json.load(file)
    os.makedirs(download_path, exist_ok=True)
    print(f"Models will be downloaded to: {download_path}")
    print("You can change the download location by using the --path argument.")
    
    for key, model_info in models.items():
        url = model_info['url']
        expected_hash = model_info['hash']
        file_path = os.path.join(download_path, os.path.basename(url))
        
        if os.path.exists(file_path) and verify_file_hash(file_path, expected_hash):
            print(f"{key} already downloaded and verified.")
            user_input = input("File already exists and is verified. Do you want to re-download it? (yes/no): ")
            if user_input.lower() != 'yes':
                continue
        
        print(f"Downloading {key} to {file_path}...")
        urllib.request.urlretrieve(url, file_path, reporthook=download_progress)
        if verify_file_hash(file_path, expected_hash):
            print(f"Successfully verified the hash for {key}.")
        else:
            raise ValueError(f"Hash mismatch for {key}, download might be corrupted.")
        os.environ[f"VIRUSNET_{key.upper()}"] = file_path

def download_progress(block_num, block_size, total_size):
    downloaded = block_num * block_size
    if total_size > 0:
        progress_percentage = min(downloaded * 100 / total_size, 100)  # Ensure percentage does not exceed 100%
        progress_bar = f"[{'=' * int(progress_percentage // 2)}{' ' * (50 - int(progress_percentage // 2))}]"
        print(f"\rDownloading: {progress_bar} {progress_percentage:.2f}%", end='')
        if downloaded >= total_size:
            print()

def verify_file_hash(file_path, expected_hash):
    sha256 = hashlib.sha256()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            sha256.update(chunk)
    calculated_hash = sha256.hexdigest()
    return calculated_hash == expected_hash

def check_files(download_path):
    models_path = os.path.join(sys.prefix, 'bin', 'models.json')
    if not os.path.exists(models_path):
        raise FileNotFoundError(f"models.json not found at {models_path}")
    
    with open(models_path, 'r') as file:
        models = json.load(file)
    for key, model_info in models.items():
        file_path = os.path.join(download_path, os.path.basename(model_info['url']))
        if not os.path.exists(file_path) or not verify_file_hash(file_path, model_info['hash']):
            return False
    return True

def run_prediction(input, output, model_paths, step_size=1000, batch_size=100, mode='binary', metagenome=False):
    """
    Function to run the R script for virus prediction using the specified arguments.
    """
    if mode == 'binary':
        if metagenome:
            r_script_path = os.path.join(os.path.dirname(__file__), "predict_binary_metagenome.r")
        else:
            r_script_path = os.path.join(os.path.dirname(__file__), "predict_binary.r")
    elif mode == 'genus':
        r_script_path = os.path.join(os.path.dirname(__file__), "predict_genus.r")
    else:
        raise ValueError(f"Invalid mode: {mode}")

    command = [
        "Rscript", r_script_path,
        '--input', input,
        '--output', output,
        '--model_binary', model_paths['binary_model'], 
        '--model_genus', model_paths['genus_model'], 
        '--labels_genus', model_paths['genus_labels'],
        '--step_size', str(step_size),
        '--batch_size', str(batch_size)
    ]
    subprocess.run(command)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='VirusNet Tool')
    subparsers = parser.add_subparsers(dest='command')

    download_parser = subparsers.add_parser('download')
    download_parser.add_argument('--path', type=str, default=os.path.expanduser('~/.virusnet'))

    predict_parser = subparsers.add_parser('predict')
    predict_parser.add_argument('--input', type=str, required=True)
    predict_parser.add_argument('--output', type=str, required=True)
    predict_parser.add_argument('--path', type=str, default=os.path.expanduser('~/.virusnet'))
    predict_parser.add_argument('--step_size', type=int, default=1000, help='Step size for prediction')
    predict_parser.add_argument('--batch_size', type=int, default=100, help='Batch size for prediction')
    predict_parser.add_argument('--mode', type=str, choices=['binary', 'genus'], default='binary', help='Prediction mode')
    predict_parser.add_argument('--metagenome', action='store_true', help='Enable metagenome mode (only applicable for binary mode)')
    
    args = parser.parse_args()

    if args.command == 'download':
        download_models(args.path)
    elif args.command == 'predict':
        if check_files(args.path):
            model_paths = {key: os.path.join(args.path, os.path.basename(info['url'])) for key, info in json.load(open(os.path.join(sys.prefix, 'bin', 'models.json'))).items()}
            run_prediction(args.input, args.output, model_paths, args.step_size, args.batch_size, args.mode, args.metagenome)
        else:
            print("Model files are missing or corrupted. Please download them again.")