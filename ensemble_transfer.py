# scripts/convert_final_edges.py (Phiên bản Nâng cấp)
import pandas as pd
import mygene
import os
import logging
from tqdm import tqdm

# --- CẤU HÌNH ---
INPUT_EDGES_PATH = 'data/edges/gene_mirna.csv'
OUTPUT_EDGES_PATH = 'data/edges/gene_mirna_ensembl.csv'
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- BẢN ĐỒ SỬA LỖI ---
# Chúng ta thêm các cặp (Tên lỗi: Ensembl ID) đã tìm được ở trên vào đây
MANUAL_CORRECTION_MAP = {
    'TMEM104': 'ENSG00000169022',
    'FAM104A': 'ENSG00000183363', # Alias for C3orf49
    'C5orf51': 'ENSG00000164419', # Alias for CFAP99
    'KIAA0895L': 'ENSG00000196733',# Alias for KATNIP
    'UHRF1BP1': 'ENSG00000148735', # Alias for ICF45
    'C19orf54': 'ENSG00000179969', # Alias for ELS1
    'NDUFA4': 'ENSG00000186518',
    'KIAA0895': 'ENSG00000137812', # Alias for CEP104
    'EFCAB2': 'ENSG00000147879'  # Alias for WDR37
}

# --- HÀM THỰC THI ---

def convert_symbols_to_ensembl_with_fallback(gene_symbols):
    """Chuyển đổi ID, sử dụng bản đồ sửa lỗi nếu mygene thất bại."""
    logging.info(f"Querying MyGene.info for {len(gene_symbols)} unique gene symbols...")
    mg = mygene.MyGeneInfo()
    
    results = mg.querymany(gene_symbols, scopes='symbol', fields='ensembl.gene', species='human', verbose=False)
    
    mapping = {}
    not_found_initially = []
    
    for res in tqdm(results, desc="Building symbol-to-ensembl map (Phase 1: Auto)"):
        query = res.get('query')
        ensembl_info = res.get('ensembl')
        if ensembl_info:
            if isinstance(ensembl_info, list):
                ensembl_id = ensembl_info[0].get('gene')
            else:
                ensembl_id = ensembl_info.get('gene')
            if ensembl_id:
                mapping[query] = ensembl_id
        else:
            not_found_initially.append(query)
    
    # Giai đoạn 2: Sửa lỗi thủ công
    if not_found_initially:
        logging.warning(f"{len(not_found_initially)} symbols not found by mygene. Attempting manual correction...")
        final_not_found = []
        for symbol in not_found_initially:
            if symbol in MANUAL_CORRECTION_MAP:
                mapping[symbol] = MANUAL_CORRECTION_MAP[symbol]
                logging.info(f"  > Manually mapped '{symbol}' to '{MANUAL_CORRECTION_MAP[symbol]}'")
            else:
                final_not_found.append(symbol)
        
        if final_not_found:
            logging.warning(f"{len(final_not_found)} symbols still could not be converted: {final_not_found}")

    logging.info(f"Successfully created mapping for {len(mapping)} genes.")
    return mapping

def convert_edge_file_to_ensembl(input_path, output_path):
    """Đọc file cạnh, chuyển đổi gene ID, và lưu file mới."""
    try:
        df = pd.read_csv(input_path)
    except FileNotFoundError:
        logging.error(f"FATAL: Input file not found at {input_path}.")
        return

    unique_symbols = df['gene_id'].unique().tolist()
    symbol_to_ensembl_map = convert_symbols_to_ensembl_with_fallback(unique_symbols)
    
    logging.info("Applying mapping to the dataframe...")
    df['gene_id'] = df['gene_id'].map(symbol_to_ensembl_map)
    
    original_rows = len(df)
    df.dropna(subset=['gene_id'], inplace=True)
    new_rows = len(df)
    
    if original_rows > new_rows:
        logging.warning(f"Dropped {original_rows - new_rows} rows due to missing Ensembl ID mapping (if any).")
        
    df.to_csv(output_path, index=False)
    logging.info(f"--- CONVERSION COMPLETE ---")
    logging.info(f"Successfully saved Ensembl-based edge file with {len(df)} rows to: {output_path}")

if __name__ == "__main__":
    convert_edge_file_to_ensembl(INPUT_EDGES_PATH, OUTPUT_EDGES_PATH)