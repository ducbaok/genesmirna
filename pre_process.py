# scripts/preprocess_interaction_data.py
import pandas as pd
import os
import logging

# --- CẤU HÌNH ---

# 1. Đường dẫn đến các file thô
RAW_DATA_DIR = 'data/raw/'
MI_RTARBASE_RAW_PATH = os.path.join(RAW_DATA_DIR, 'hsa_MTI_homo.csv')
TARGETSCAN_RAW_PATH = os.path.join(RAW_DATA_DIR, 'Predicted_Targets_Info.txt')
# File mapping mới tải về
FAMILY_INFO_PATH = os.path.join(RAW_DATA_DIR, 'miR_Family_Info.txt') 

# 2. Đường dẫn output
PROCESSED_DATA_DIR = 'data/processed/'
MI_RTARBASE_PROCESSED_PATH = os.path.join(PROCESSED_DATA_DIR, 'mirtarbase_processed.csv')
TARGETSCAN_PROCESSED_PATH = os.path.join(PROCESSED_DATA_DIR, 'targetscan_processed.csv')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- CÁC HÀM XỬ LÝ ---

def preprocess_mirtarbase(raw_path, processed_path):
    """Tiền xử lý miRTarBase (Giữ nguyên logic cũ)."""
    logging.info(f"Processing miRTarBase file from: {raw_path}")
    try:
        df = pd.read_csv(raw_path)
        
        # 1. Lọc loài người & Chọn cột
        df_human = df[df['Species (miRNA)'] == 'hsa']
        df_processed = df_human[['miRNA', 'Target Gene']].copy()
        df_processed.rename(columns={'miRNA': 'mirna_id', 'Target Gene': 'gene_id'}, inplace=True)
        
        # 2. Làm sạch
        df_processed.drop_duplicates(inplace=True)
        
        os.makedirs(os.path.dirname(processed_path), exist_ok=True)
        df_processed.to_csv(processed_path, index=False)
        logging.info(f"Saved miRTarBase data to: {processed_path}")

    except FileNotFoundError:
        logging.error(f"FATAL: miRTarBase file not found at {raw_path}.")
    except Exception as e:
        logging.error(f"Error processing miRTarBase: {e}")

def load_family_mapping(family_info_path): # <--- HÀM MỚI QUAN TRỌNG
    """
    Tạo từ điển ánh xạ từ 'miR Family' sang danh sách các 'MiRBase ID' cụ thể (có hsa-).
    """
    logging.info(f"Loading family mapping from: {family_info_path}")
    try:
        # Đọc file txt tab-separated
        df = pd.read_csv(family_info_path, sep='\t')
        
        # 1. Lọc loài người (Species ID = 9606)
        # Cột 'Species ID' trong file là số nguyên
        df_human = df[df['Species ID'] == 9606]
        
        # 2. Tạo Map: Family -> List of IDs (hsa-miR-...)
        # Cột 'miR family' là key (VD: let-7/98), 'MiRBase ID' là value (VD: hsa-let-7a-5p)
        mapping = df_human.groupby('miR family')['MiRBase ID'].apply(list).to_dict()
        
        logging.info(f"Loaded mapping for {len(mapping)} human miRNA families.")
        return mapping
    except Exception as e:
        logging.error(f"FATAL: Could not load family info. Check path and format. Error: {e}")
        return {}

def preprocess_targetscan(raw_path, processed_path, family_mapping):
    """
    Tiền xử lý TargetScan với logic 'Bùng nổ' (Explode) sử dụng map.
    """
    logging.info(f"Processing TargetScan file from: {raw_path}")
    try:
        # Đọc file raw (thường là file tab-separated)
        df = pd.read_csv(raw_path, sep='\t')
        
        # Kiểm tra tên cột (TargetScan thường dùng 'miR Family' và 'Gene Symbol')
        if 'miR Family' not in df.columns or 'Gene Symbol' not in df.columns:
            logging.error("FATAL: TargetScan columns mismatch. Expected 'miR Family' and 'Gene Symbol'.")
            return

        expanded_rows = []
        unique_pairs = set() # Để loại bỏ trùng lặp nhanh
        
        logging.info("Exploding miRNA families... (This may take a minute)")
        
        # Duyệt qua từng dòng dự đoán
        count_mapped = 0
        count_missed = 0
        
        for _, row in df.iterrows():
            family_key = row['miR Family']
            gene = row['Gene Symbol']
            
            # Tra cứu trong từ điển mapping
            if family_key in family_mapping:
                specific_mirnas = family_mapping[family_key]
                for mirna_id in specific_mirnas:
                    # mirna_id lúc này đã chuẩn dạng: 'hsa-miR-122-5p'
                    pair = (mirna_id, gene)
                    if pair not in unique_pairs:
                        expanded_rows.append({'mirna_id': mirna_id, 'gene_id': gene})
                        unique_pairs.add(pair)
                count_mapped += 1
            else:
                count_missed += 1
                # Tùy chọn: Bạn có thể log các family không tìm thấy map nếu cần
                # logging.debug(f"Missing map for family: {family_key}")

        # Tạo DataFrame kết quả
        df_processed = pd.DataFrame(expanded_rows)
        
        # Lưu file
        os.makedirs(os.path.dirname(processed_path), exist_ok=True)
        df_processed.to_csv(processed_path, index=False)
        
        logging.info(f"TargetScan Processing Complete.")
        logging.info(f"  - Original predicted families processed: {count_mapped}")
        logging.info(f"  - Families skipped (no human map): {count_missed}")
        logging.info(f"  - Final expanded interactions: {len(df_processed)}")
        logging.info(f"Saved to: {processed_path}")

    except FileNotFoundError:
        logging.error(f"FATAL: Raw TargetScan file not found at {raw_path}.")
    except Exception as e:
        logging.error(f"An error occurred while processing TargetScan: {e}")

if __name__ == "__main__":
    logging.info("--- Starting Enhanced Preprocessing Pipeline ---")
    
    os.makedirs(RAW_DATA_DIR, exist_ok=True)
    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
    
    # 1. miRTarBase
    preprocess_mirtarbase(MI_RTARBASE_RAW_PATH, MI_RTARBASE_PROCESSED_PATH)
    
    # 2. TargetScan (Cần load map trước)
    if os.path.exists(FAMILY_INFO_PATH):
        fam_map = load_family_mapping(FAMILY_INFO_PATH)
        if fam_map:
            preprocess_targetscan(TARGETSCAN_RAW_PATH, TARGETSCAN_PROCESSED_PATH, fam_map)
    else:
        logging.error(f"FATAL: Family Info file missing at {FAMILY_INFO_PATH}. Cannot process TargetScan correctly.")
    
    logging.info("--- Pipeline Finished ---")