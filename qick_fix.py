# scripts/quick_fix_mirna.py
import pandas as pd
import os
import logging

# --- CẤU HÌNH ---
INPUT_PATH = 'data/features/mirnas.tsv'  # File cũ của bạn
OUTPUT_PATH = 'data/features/mirnas.tsv' # Ghi đè lại chính nó (hoặc đổi tên nếu muốn backup)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

def normalize_mirna_id(mir_id):
    """
    Chuẩn hóa ID: 
    1. hsa-let-7a-1 -> hsa-let-7a (Bỏ hậu tố copy number)
    2. hsa-mir-122 -> hsa-miR-122 (Chuẩn hóa chữ hoa/thường)
    """
    # Chuyển 'mir' thành 'miR' (theo chuẩn TargetScan/miRBase)
    new_id = mir_id.replace('mir', 'miR')
    
    # Bỏ hậu tố số bản sao (-1, -2) nếu có
    # VD: hsa-let-7a-1 -> hsa-let-7a
    parts = new_id.split('-')
    if parts[-1].isdigit():
        new_id = '-'.join(parts[:-1])
        
    return new_id

def fix_existing_matrix():
    logging.info(f"Reading existing matrix from {INPUT_PATH}...")
    try:
        # Load dữ liệu cũ
        df = pd.read_csv(INPUT_PATH, sep='\t', index_col=0)
        original_count = len(df)
        
        logging.info("Normalizing IDs...")
        # Áp dụng hàm chuẩn hóa cho Index
        df.index = df.index.map(normalize_mirna_id)
        
        # Gom nhóm các dòng bị trùng tên sau khi chuẩn hóa (Lấy trung bình)
        # VD: let-7a-1 và let-7a-2 giờ cùng tên let-7a -> Gộp lại
        df = df.groupby(df.index).mean()
        
        new_count = len(df)
        logging.info(f"Done. Reduced from {original_count} to {new_count} unique miRNAs.")
        
        # Lưu lại
        df.to_csv(OUTPUT_PATH, sep='\t')
        logging.info(f"Saved fixed matrix to {OUTPUT_PATH}")
        
    except FileNotFoundError:
        logging.error(f"FATAL: Could not find {INPUT_PATH}. Please re-download data from GDC if you lost this file.")

if __name__ == "__main__":
    fix_existing_matrix()