# scripts/pre_process.py
import pandas as pd
import os
import logging

# --- CẤU HÌNH ---
RAW_DATA_DIR = 'data/raw/'
MI_RTARBASE_RAW_PATH = os.path.join(RAW_DATA_DIR, 'hsa_MTI_homo.csv')
TARGETSCAN_RAW_PATH = os.path.join(RAW_DATA_DIR, 'Predicted_Targets_Info.txt')
FAMILY_INFO_PATH = os.path.join(RAW_DATA_DIR, 'miR_Family_Info.txt') # File mapping mới

PROCESSED_DATA_DIR = 'data/processed/'
MI_RTARBASE_PROCESSED_PATH = os.path.join(PROCESSED_DATA_DIR, 'mirtarbase_processed.csv')
TARGETSCAN_PROCESSED_PATH = os.path.join(PROCESSED_DATA_DIR, 'targetscan_processed.csv')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def preprocess_mirtarbase(raw_path, processed_path):
    logging.info(f"Processing miRTarBase: {raw_path}")
    try:
        df = pd.read_csv(raw_path)
        df_human = df[df['Species (miRNA)'] == 'hsa']
        df_processed = df_human[['miRNA', 'Target Gene']].copy()
        df_processed.rename(columns={'miRNA': 'mirna_id', 'Target Gene': 'gene_id'}, inplace=True)
        df_processed.drop_duplicates(inplace=True)
        os.makedirs(os.path.dirname(processed_path), exist_ok=True)
        df_processed.to_csv(processed_path, index=False)
        logging.info(f"Saved miRTarBase: {len(df_processed)} interactions.")
    except Exception as e:
        logging.error(f"Error miRTarBase: {e}")

def load_family_mapping(family_info_path):
    """Map: Family -> List of Specific IDs (hsa-miR-...)"""
    logging.info(f"Loading family map: {family_info_path}")
    try:
        df = pd.read_csv(family_info_path, sep='\t')
        df_human = df[df['Species ID'] == 9606] # Human only
        # Group by Family -> List of MiRBase IDs
        mapping = df_human.groupby('miR family')['MiRBase ID'].apply(list).to_dict()
        return mapping
    except Exception as e:
        logging.error(f"Error loading family info: {e}")
        return {}

def preprocess_targetscan(raw_path, processed_path, family_mapping):
    logging.info(f"Processing TargetScan with Family Expansion: {raw_path}")
    try:
        df = pd.read_csv(raw_path, sep='\t')
        expanded_rows = []
        unique_pairs = set()
        
        for _, row in df.iterrows():
            family = row['miR Family']
            gene = row['Gene Symbol']
            
            # Tra cứu danh sách con trong họ
            specific_mirnas = family_mapping.get(family, [])
            
            # Nếu không thấy trong map (hiếm), giữ nguyên tên family (có thể sẽ bị lọc sau)
            if not specific_mirnas:
                # Cố gắng sửa tên cơ bản nếu thiếu hsa-
                if not family.startswith('hsa-') and 'miR' in family:
                    specific_mirnas = ['hsa-' + family]
                else:
                    specific_mirnas = [family]

            for mirna_id in specific_mirnas:
                pair = (mirna_id, gene)
                if pair not in unique_pairs:
                    expanded_rows.append({'mirna_id': mirna_id, 'gene_id': gene})
                    unique_pairs.add(pair)

        df_final = pd.DataFrame(expanded_rows)
        os.makedirs(os.path.dirname(processed_path), exist_ok=True)
        df_final.to_csv(processed_path, index=False)
        logging.info(f"Saved Expanded TargetScan: {len(df_final)} interactions.")
    except Exception as e:
        logging.error(f"Error TargetScan: {e}")

if __name__ == "__main__":
    os.makedirs(RAW_DATA_DIR, exist_ok=True)
    os.makedirs(PROCESSED_DATA_DIR, exist_ok=True)
    
    preprocess_mirtarbase(MI_RTARBASE_RAW_PATH, MI_RTARBASE_PROCESSED_PATH)
    
    if os.path.exists(FAMILY_INFO_PATH):
        fam_map = load_family_mapping(FAMILY_INFO_PATH)
        preprocess_targetscan(TARGETSCAN_RAW_PATH, TARGETSCAN_PROCESSED_PATH, fam_map)
    else:
        logging.error("Missing miR_Family_Info.txt. Cannot expand TargetScan families properly.")