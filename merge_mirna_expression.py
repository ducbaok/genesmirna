# scripts/merge_mirna_expression.py
import pandas as pd
import os
from tqdm import tqdm
import json
import logging

# --- CẤU HÌNH ---
DOWNLOADED_FILES_DIR = 'miRNA_expression' # Thư mục chứa file GDC tải về
MANIFEST_FILE_PATH = os.path.join(DOWNLOADED_FILES_DIR, 'MANIFEST.txt')
METADATA_FILE_PATH = os.path.join(DOWNLOADED_FILES_DIR, 'METADATA.json')
OUTPUT_MATRIX_PATH = 'data/features/mirnas.tsv'

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def create_file_to_patient_map(metadata_path):
    logging.info(f"Parsing metadata from: {metadata_path}")
    mapping = {}
    try:
        with open(metadata_path, 'r') as f:
            metadata = json.load(f)
        for item in metadata:
            file_name = item.get('file_name')
            if item.get('associated_entities'):
                submitter_id = item['associated_entities'][0].get('entity_submitter_id')
                if file_name and submitter_id:
                    patient_id = '-'.join(submitter_id.split('-')[:3])
                    mapping[file_name] = patient_id
        return mapping
    except Exception as e:
        logging.error(f"Error loading metadata: {e}")
        return None

def normalize_mirna_id(mir_id):
    """Chuẩn hóa ID: hsa-let-7a-1 -> hsa-let-7a, mir -> miR"""
    # 1. Chuyển 'mir' thành 'miR' (chuẩn chung)
    new_id = mir_id.replace('mir', 'miR')
    # 2. Bỏ hậu tố số bản sao (-1, -2)
    parts = new_id.split('-')
    if parts[-1].isdigit():
        new_id = '-'.join(parts[:-1])
    return new_id

def merge_mirna_files(data_dir, manifest_path, metadata_path, output_path):
    logging.info("Starting miRNA file merge process...")
    file_to_patient_map = create_file_to_patient_map(metadata_path)
    if not file_to_patient_map: return

    try:
        manifest_df = pd.read_csv(manifest_path, sep='\t')
    except Exception as e:
        logging.error(f"Manifest error: {e}")
        return

    all_mirna_series = []
    for _, row in tqdm(manifest_df.iterrows(), total=manifest_df.shape[0], desc="Processing files"):
        file_name = row['filename']
        file_path = os.path.join(data_dir, file_name)
        
        if not os.path.exists(file_path): continue
        
        patient_id = file_to_patient_map.get(os.path.basename(file_name))
        if not patient_id: continue

        try:
            df = pd.read_csv(file_path, sep='\t', usecols=['miRNA_ID', 'reads_per_million_miRNA_mapped'], index_col='miRNA_ID')
            series = df['reads_per_million_miRNA_mapped']
            series.name = patient_id
            all_mirna_series.append(series)
        except Exception:
            pass

    if not all_mirna_series:
        logging.error("No data processed.")
        return

    logging.info("Merging matrix...")
    final_matrix = pd.concat(all_mirna_series, axis=1).fillna(0)
    
    # Gom nhóm các bệnh nhân trùng lặp (nếu có)
    final_matrix = final_matrix.groupby(final_matrix.columns, axis=1).mean()
    
    # --- BƯỚC CHUẨN HÓA QUAN TRỌNG ---
    logging.info("Normalizing miRNA IDs (Stem -> Base)...")
    final_matrix.index = final_matrix.index.map(normalize_mirna_id)
    # Gom nhóm các dòng trùng tên sau khi chuẩn hóa (VD: let-7a-1 và let-7a-2 -> let-7a)
    final_matrix = final_matrix.groupby(final_matrix.index).mean()
    
    # Lọc các dòng toàn 0
    final_matrix = final_matrix.loc[(final_matrix.sum(axis=1) > 0)]

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    final_matrix.to_csv(output_path, sep='\t', index_label='miRNA_ID')
    logging.info(f"Saved merged data to {output_path}. Shape: {final_matrix.shape}")

if __name__ == "__main__":
    merge_mirna_files(DOWNLOADED_FILES_DIR, MANIFEST_FILE_PATH, METADATA_FILE_PATH, OUTPUT_MATRIX_PATH)