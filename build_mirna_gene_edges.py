import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm
import os
import logging
import re

# --- CONFIG ---
GENE_PATH = 'data/features/genes_expr.txt'
MIR_PATH = 'data/features/mirnas.tsv'
MIRTAR_PATH = 'data/processed/mirtarbase_processed.csv'
TARGET_PATH = 'data/processed/targetscan_processed.csv'
OUTPUT_PATH = 'data/edges/gene_mirna.csv'

P_THRESH = 0.05
R_THRESH = -0.1 # Tương quan nghịch
BONUS = 0.1

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')

def clean_cols(df):
    df.columns = ['-'.join(x.split('-')[:3]) for x in df.columns]
    return df.groupby(df.columns, axis=1).mean()

def build_edges():
    logging.info("Loading data...")
    try:
        # Load và Clean Data
        gene_df = pd.read_csv(GENE_PATH, sep='\t', index_col=0)
        if 'Entrez_Gene_Id' in gene_df.columns: gene_df.drop(columns=['Entrez_Gene_Id'], inplace=True)
        gene_df = clean_cols(gene_df[~gene_df.index.duplicated()])
        
        mir_df = pd.read_csv(MIR_PATH, sep='\t', index_col=0)
        mir_df = clean_cols(mir_df)
        
        # Intersection Patients
        common = gene_df.columns.intersection(mir_df.columns)
        if len(common) == 0: return logging.error("No common patients.")
        gene_df, mir_df = gene_df[common], mir_df[common]
        
        # Load Candidates
        val = pd.read_csv(MIRTAR_PATH); val['validated'] = True
        pred = pd.read_csv(TARGET_PATH); pred['validated'] = False
        candidates = pd.concat([val, pred]).drop_duplicates(subset=['mirna_id', 'gene_id'])
        
        # --- BƯỚC 1: TẠO MAP (Regex) ---
        # Map: Core Name (mir-122) -> Precursor ID (hsa-mir-122)
        expr_map = {}
        for idx in mir_df.index:
            # Regex bắt: (mir hoặc let) - (chuỗi số/chữ)
            m = re.search(r'(mir|let)-([0-9a-z]+)', idx, re.IGNORECASE)
            if m:
                # Key chuẩn hóa về dạng: mir-122 (chữ thường)
                core = m.group(0).lower().replace('mir', 'mir')
                expr_map[core] = idx # Lưu ID gốc trong file
                
        logging.info(f"Expression Map created for {len(expr_map)} precursors.")

        # --- BƯỚC 2: TÍNH TOÁN ---
        results = []
        gene_cache = gene_df.to_dict('index')
        mir_cache = mir_df.to_dict('index')
        valid_genes = set(gene_df.index)

        for _, row in tqdm(candidates.iterrows(), total=len(candidates), desc="Processing"):
            mature_id = row['mirna_id'] # VD: hsa-miR-122-5p
            gene = row['gene_id']
            
            if gene not in valid_genes: continue

            # Tìm ID Precursor tương ứng (Cha)
            target_pre_id = None
            m = re.search(r'(mir|let)-([0-9a-z]+)', mature_id, re.IGNORECASE)
            if m:
                core = m.group(0).lower().replace('mir', 'mir')
                target_pre_id = expr_map.get(core)
            
            if not target_pre_id: continue # Không có data biểu hiện -> Bỏ qua

            # Lấy vector & Tính Correlation
            try:
                vec_m = pd.Series(mir_cache[target_pre_id])
                vec_g = pd.Series(gene_cache[gene])
                
                if vec_m.var() == 0 or vec_g.var() == 0: continue
                
                r, p = pearsonr(vec_m, vec_g)
                
                if p < P_THRESH and r < R_THRESH:
                    w = abs(r)
                    if row['validated']: w = min(1.0, w + BONUS)
                    
                    # --- QUAN TRỌNG: LƯU ID PRECURSOR ---
                    results.append({
                        'mirna_id': target_pre_id, # Lưu ID cha để khớp với Node Features
                        'gene_id': gene,
                        'weight': w
                    })
            except: continue

        # --- BƯỚC 3: GỘP TRÙNG LẶP ---
        # Nếu cả 5p và 3p cùng trỏ vào 1 gene -> Giữ cái có trọng số cao nhất
        final_df = pd.DataFrame(results)
        if not final_df.empty:
            final_df = final_df.sort_values('weight', ascending=False)
            final_df = final_df.drop_duplicates(subset=['mirna_id', 'gene_id'], keep='first')
            
            os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
            final_df.to_csv(OUTPUT_PATH, index=False)
            logging.info(f"SUCCESS: {len(final_df)} edges saved to {OUTPUT_PATH}")
        else:
            logging.warning("No edges passed filter.")

    except Exception as e:
        logging.error(f"Pipeline failed: {e}")

if __name__ == "__main__":
    build_edges()