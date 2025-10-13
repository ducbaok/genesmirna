# scripts/preprocess_interaction_data.py
import pandas as pd
import os
import logging

# --- CẤU HÌNH ---

# 1. Đường dẫn đến các file thô bạn đã tải về
#    (Giả sử bạn đã đặt chúng vào một thư mục 'data/raw/')
RAW_DATA_DIR = 'data/raw/'
MI_RTARBASE_RAW_PATH = os.path.join(RAW_DATA_DIR, 'hsa_MTI_homo.csv')
TARGETSCAN_RAW_PATH = os.path.join(RAW_DATA_DIR, 'Predicted_Targets_Info.txt')

# 2. Đường dẫn đến thư mục chứa các file đã được xử lý
PROCESSED_DATA_DIR = 'data/processed/'
MI_RTARBASE_PROCESSED_PATH = os.path.join(PROCESSED_DATA_DIR, 'mirtarbase_processed.csv')
TARGETSCAN_PROCESSED_PATH = os.path.join(PROCESSED_DATA_DIR, 'targetscan_processed.csv')

# Cấu hình logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- HÀM THỰC THI ---

def preprocess_mirtarbase(raw_path, processed_path):
    """
    Tiền xử lý file miRTarBase từ định dạng CSV.
    """
    logging.info(f"Processing miRTarBase file from: {raw_path}")
    try:
        # --- THAY ĐỔI QUAN TRỌNG Ở ĐÂY ---
        # Đọc file CSV thay vì Excel
        df = pd.read_csv(raw_path)
        # ---------------------------------
        
        # 1. Lọc loài người
        df_human = df[df['Species (miRNA)'] == 'hsa']
        
        # 2. Chọn và đổi tên cột
        df_processed = df_human[['miRNA', 'Target Gene']].copy()
        df_processed.rename(columns={
            'miRNA': 'mirna_id',
            'Target Gene': 'gene_id'
        }, inplace=True)
        
        # 3. Loại bỏ các dòng trùng lặp
        df_processed.drop_duplicates(inplace=True)
        
        # 4. Lưu file
        os.makedirs(os.path.dirname(processed_path), exist_ok=True)
        df_processed.to_csv(processed_path, index=False)
        
        logging.info(f"Successfully processed {len(df_processed)} unique human interactions from miRTarBase.")
        logging.info(f"Saved to: {processed_path}")

    except FileNotFoundError:
        logging.error(f"FATAL: Raw miRTarBase file not found at {raw_path}. Please check the file name and path.")
    except Exception as e:
        logging.error(f"An error occurred while processing miRTarBase: {e}")


def preprocess_targetscan(raw_path, processed_path):
    """
    Tiền xử lý file TargetScan:
    1. Tải file .txt (dạng tab-separated).
    2. Chọn và đổi tên các cột cần thiết.
    3. Lưu dưới dạng file CSV sạch.
    """
    logging.info(f"Processing TargetScan file from: {raw_path}")
    try:
        df = pd.read_csv(raw_path, sep='\t')
        
        # 1. Chọn và đổi tên cột
        # Chúng ta dùng 'Gene Symbol' để nhất quán với file genes_expr.tsv
        df_processed = df[['miR Family', 'Gene Symbol']].copy()
        df_processed.rename(columns={
            'miR Family': 'mirna_id',
            'Gene Symbol': 'gene_id'
        }, inplace=True)

        # 2. Loại bỏ các dòng trùng lặp
        df_processed.drop_duplicates(inplace=True)
        
        # 3. Lưu file
        os.makedirs(os.path.dirname(processed_path), exist_ok=True)
        df_processed.to_csv(processed_path, index=False)
        
        logging.info(f"Successfully processed {len(df_processed)} unique predicted interactions.")
        logging.info(f"Saved to: {processed_path}")

    except FileNotFoundError:
        logging.error(f"FATAL: Raw TargetScan file not found at {raw_path}. Please download and unzip it first.")
    except Exception as e:
        logging.error(f"An error occurred while processing TargetScan: {e}")


if __name__ == "__main__":
    logging.info("--- Starting Preprocessing Pipeline for Interaction Data ---")
    
    # Tạo các thư mục cần thiết
    os.makedirs(RAW_DATA_DIR, exist_ok=True)
    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
    
    # Chạy tiền xử lý cho từng file
    preprocess_mirtarbase(MI_RTARBASE_RAW_PATH, MI_RTARBASE_PROCESSED_PATH)
    preprocess_targetscan(TARGETSCAN_RAW_PATH, TARGETSCAN_PROCESSED_PATH)
    
    logging.info("--- Preprocessing Pipeline Finished ---")