# scripts/merge_mirna_files.py
import pandas as pd
import os
from tqdm import tqdm
import json
import logging

# --- CẤU HÌNH ---
# 1. Đường dẫn đến thư mục chứa dữ liệu miRNA đã tải về từ GDC
DOWNLOADED_FILES_DIR = 'miRNA_expression'
# 2. Đường dẫn đến file manifest bạn đã tải từ GDC
MANIFEST_FILE_PATH = os.path.join(DOWNLOADED_FILES_DIR, 'MANIFEST.txt') # Cập nhật tên file manifest cho đúng
# 3. Đường dẫn đến file metadata.json
METADATA_FILE_PATH = os.path.join(DOWNLOADED_FILES_DIR, 'METADATA.json') # Cập nhật tên file metadata cho đúng
# 4. Đường dẫn file output cuối cùng
OUTPUT_MATRIX_PATH = 'data/features/mirnas.tsv'



# Cấu hình logging để dễ dàng gỡ lỗi
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- HÀM THỰC THI ---

def create_file_to_patient_map(metadata_path):
    """Tạo dictionary ánh xạ từ tên file (filename) sang mã bệnh nhân đã làm sạch."""
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
        
        logging.info(f"Successfully created map for {len(mapping)} files.")
        return mapping
    except FileNotFoundError:
        logging.error(f"FATAL: Metadata file not found at {metadata_path}. Please check the path.")
        return None
    except json.JSONDecodeError:
        logging.error(f"FATAL: Metadata file {metadata_path} is not a valid JSON file.")
        return None

def merge_mirna_files(data_dir, manifest_path, metadata_path, output_path):
    """Gộp nhiều file định lượng miRNA từ GDC thành một ma trận duy nhất."""
    logging.info("Starting miRNA file merge process...")

    file_to_patient_map = create_file_to_patient_map(metadata_path)
    if not file_to_patient_map:
        return

    try:
        manifest_df = pd.read_csv(manifest_path, sep='\t')
    except FileNotFoundError:
        logging.error(f"FATAL: Manifest file not found at {manifest_path}. Please check the path.")
        return

    all_mirna_series = []
    
    for _, row in tqdm(manifest_df.iterrows(), total=manifest_df.shape[0], desc="Processing GDC files"):
        # file_id = row['id'] # Không cần dùng file_id để xây dựng đường dẫn nữa
        file_name_from_manifest = row['filename']
        
        # --- THAY ĐỔI QUAN TRỌNG Ở ĐÂY ---
        # Xây dựng đường dẫn chỉ từ thư mục cha và tên file trong manifest
        # vì tên file đã chứa sẵn đường dẫn tương đối.
        file_path = os.path.join(data_dir, file_name_from_manifest)
        # ------------------------------------

        if not os.path.exists(file_path):
            logging.warning(f"File not found at {file_path}. Skipping.")
            continue
            
        # Lấy tên file thực tế (không có đường dẫn) để tra cứu mã bệnh nhân
        actual_file_name = os.path.basename(file_path)
        patient_id = file_to_patient_map.get(actual_file_name)
        if not patient_id:
            logging.warning(f"Could not find patient ID for file {actual_file_name}. Skipping.")
            continue

        try:
            df = pd.read_csv(file_path, sep='\t', usecols=['miRNA_ID', 'reads_per_million_miRNA_mapped'], 
                             index_col='miRNA_ID')
            series = df['reads_per_million_miRNA_mapped']
            series.name = patient_id
            all_mirna_series.append(series)

        except Exception as e:
            logging.error(f"Could not process file {file_path}: {e}")

    if not all_mirna_series:
        logging.error("No data was processed. Please check your file paths and formats.")
        return

    logging.info(f"Merging data from {len(all_mirna_series)} files into a single matrix...")
    final_matrix = pd.concat(all_mirna_series, axis=1)
    final_matrix = final_matrix.fillna(0)
    final_matrix = final_matrix.groupby(final_matrix.columns, axis=1).mean()
    final_matrix = final_matrix.loc[(final_matrix.sum(axis=1) > 0)]

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    final_matrix.to_csv(output_path, sep='\t', index_label='miRNA_ID')

    logging.info(f"--- PIPELINE FINISHED ---")
    logging.info(f"Successfully merged {len(final_matrix.columns)} patients into {output_path}")
    logging.info(f"Final matrix dimensions: {final_matrix.shape[0]} miRNAs x {final_matrix.shape[1]} patients.")


if __name__ == "__main__":
    # Cập nhật các đường dẫn dưới đây cho đúng với máy của bạn
    # DOWNLOADED_FILES_DIR = 'C:/Users/YourUser/Downloads/GDC_miRNA_Data'
    # MANIFEST_FILE_PATH = os.path.join(DOWNLOADED_FILES_DIR, 'gdc_manifest.2025-10-13.txt')
    # METADATA_FILE_PATH = os.path.join(DOWNLOADED_FILES_DIR, 'metadata.cart.2025-10-13.json')
    
    merge_mirna_files(DOWNLOADED_FILES_DIR, MANIFEST_FILE_PATH, METADATA_FILE_PATH, OUTPUT_MATRIX_PATH)