# scripts/ensemble_transfer.py
import pandas as pd
import mygene
import os
import logging
from tqdm import tqdm

# --- CẤU HÌNH ---
INPUT_EDGES_PATH = 'data/edges/gene_mirna.csv'
OUTPUT_EDGES_PATH = 'data/edges/gene_mirna_ensembl.csv'
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- BẢN ĐỒ SỬA LỖI THỦ CÔNG (FINAL VERSION) ---
# Đã bổ sung đầy đủ 18 gene cuối cùng và các gene phổ biến khác
MANUAL_CORRECTION_MAP = {
    # Các gene từ danh sách 18 gene cuối cùng của bạn
    'TRNP1': 'ENSG00000253368',
    'AGPS': 'ENSG00000018510',
    'EXT2': 'ENSG00000151348',
    'NBL1': 'ENSG00000158747',
    'KIAA1045': 'ENSG00000122733', # Alias của PHF24
    'NPC1': 'ENSG00000141458',
    'FMR1': 'ENSG00000102081',
    'PLSCR3': 'ENSG00000187838',
    'ELP4': 'ENSG00000109911',
    'EGLN2': 'ENSG00000269858',
    'PDE11A': 'ENSG00000128655',
    'ZFP91': 'ENSG00000186660',
    'ASB3': 'ENSG00000115239',
    'MEMO1': 'ENSG00000162959',
    'ZNF788': 'ENSG00000214189', # Alias của ZNF788P
    'C15orf37': 'ENSG00000259642', # Alias của ST20-AS1
    'LUZP6': 'ENSG00000267697',
    'PAK6': 'ENSG00000137843',

    # Các gene phổ biến khác (giữ lại để an toàn)
    'PHB': 'ENSG00000167085',      # PHB1
    'MARCH5': 'ENSG00000198060',   # MARCHF5
    'DDX58': 'ENSG00000107201',    # RIGI
    'KIAA0101': 'ENSG00000166803', # PCLAF
    'FAM126B': 'ENSG00000155744',  # HYCC2
    'SEPT2': 'ENSG00000168385',    # SEPTIN2
    'SEPT7': 'ENSG00000163331',    # SEPTIN7
    'H2AFX': 'ENSG00000188486',    # H2AX
    'TMEM104': 'ENSG00000169022',
    'FAM104A': 'ENSG00000183363',
    'C5orf51': 'ENSG00000164419',
    'KIAA0895L': 'ENSG00000196733',
    'UHRF1BP1': 'ENSG00000148735',
    'C19orf54': 'ENSG00000179969',
    'NDUFA4': 'ENSG00000186518',
    'KIAA0895': 'ENSG00000137812',
    'EFCAB2': 'ENSG00000147879'
}

# --- HÀM THỰC THI ---

def convert_symbols_to_ensembl_with_fallback(gene_symbols):
    """Chuyển đổi ID với phạm vi tìm kiếm mở rộng và Map thủ công."""
    logging.info(f"Querying MyGene.info for {len(gene_symbols)} unique gene symbols...")
    mg = mygene.MyGeneInfo()
    
    # Tìm kiếm mở rộng cả alias và tên cũ
    results = mg.querymany(
        gene_symbols, 
        scopes='symbol,alias,prev_symbol,reporter', 
        fields='ensembl.gene', 
        species='human', 
        verbose=False
    )
    
    mapping = {}
    not_found_initially = []
    
    for res in tqdm(results, desc="Building map"):
        query = res.get('query')
        ensembl_info = res.get('ensembl')
        
        if ensembl_info:
            if isinstance(ensembl_info, list):
                ensembl_id = ensembl_info[0].get('gene')
            else:
                ensembl_id = ensembl_info.get('gene')
            
            if ensembl_id:
                mapping[query] = ensembl_id
        
        if query not in mapping:
            not_found_initially.append(query)
            
    # Giai đoạn 2: Sửa lỗi thủ công
    final_not_found = []
    if not_found_initially:
        not_found_initially = list(set(not_found_initially))
        logging.info(f"Checking manual map for {len(not_found_initially)} missing genes...")
        
        for symbol in not_found_initially:
            if symbol in MANUAL_CORRECTION_MAP:
                mapping[symbol] = MANUAL_CORRECTION_MAP[symbol]
            else:
                final_not_found.append(symbol)
        
        if final_not_found:
            logging.warning(f"STILL MISSING {len(final_not_found)} GENES: {final_not_found}")
        else:
            logging.info("All missing genes resolved manually!")

    return mapping

def convert_edge_file_to_ensembl(input_path, output_path):
    try:
        df = pd.read_csv(input_path)
    except FileNotFoundError:
        logging.error(f"FATAL: Input file not found at {input_path}.")
        return

    unique_symbols = df['gene_id'].unique().tolist()
    symbol_to_ensembl_map = convert_symbols_to_ensembl_with_fallback(unique_symbols)
    
    logging.info("Applying mapping...")
    df['gene_id'] = df['gene_id'].map(symbol_to_ensembl_map)
    
    original_rows = len(df)
    df.dropna(subset=['gene_id'], inplace=True)
    new_rows = len(df)
    
    logging.info(f"--- DONE ---")
    logging.info(f"Saved to: {output_path}")
    logging.info(f"Rows retained: {new_rows}/{original_rows} ({new_rows/original_rows*100:.1f}%)")
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)

if __name__ == "__main__":
    convert_edge_file_to_ensembl(INPUT_EDGES_PATH, OUTPUT_EDGES_PATH)