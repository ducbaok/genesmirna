# scripts/build_mirna_gene_edges.py
import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm
import os
import logging

# --- CẤU HÌNH ---
GENE_EXPR_PATH = 'data/features/genes_expr.txt' 
MIRNA_EXPR_PATH = 'data/features/mirnas.tsv'

MI_RTARBASE_PATH = 'data/processed/mirtarbase_processed.csv'
TARGETSCAN_PATH = 'data/processed/targetscan_processed.csv'

OUTPUT_PATH = 'data/edges/gene_mirna.csv'
P_VALUE_THRESHOLD = 0.05
VALIDATED_BONUS = 0.1

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- HÀM THỰC THI ---

def clean_patient_ids(df):
    """Chuẩn hóa các cột tên bệnh nhân về dạng 12 ký tự (TCGA-XX-XXXX)."""
    # Tạo một hàm lambda để áp dụng cho từng tên cột
    # Tách chuỗi theo dấu '-', lấy 3 phần đầu, rồi ghép lại
    # Ví dụ: 'TCGA-05-4244-01' -> 'TCGA-05-4244'
    rename_function = lambda col: '-'.join(col.split('-')[:3])
    df = df.rename(columns=rename_function)
    return df

def load_expression_data(gene_path, mirna_path):
    """Tải và chuẩn bị dữ liệu biểu hiện với mã bệnh nhân đã được chuẩn hóa."""
    logging.info("Loading expression data...")
    try:
        gene_expr = pd.read_csv(gene_path, sep='\t', index_col='Hugo_Symbol')
        if 'Entrez_Gene_Id' in gene_expr.columns:
            gene_expr = gene_expr.drop(columns=['Entrez_Gene_Id'])
            
        mirna_expr = pd.read_csv(mirna_path, sep='\t', index_col=0)
        
        # --- THAY ĐỔI QUAN TRỌNG Ở ĐÂY ---
        # Chuẩn hóa tên cột (mã bệnh nhân) cho cả hai dataframe
        logging.info("Normalizing patient IDs...")
        gene_expr = clean_patient_ids(gene_expr)
        mirna_expr = clean_patient_ids(mirna_expr)
        # ---------------------------------
        
        common_patients = gene_expr.columns.intersection(mirna_expr.columns)
        if len(common_patients) == 0:
            logging.warning("Found 0 common patients. Please check patient ID formats in source files.")
            # Trả về DataFrame rỗng để tránh lỗi sau này, nhưng vẫn in cảnh báo
            return pd.DataFrame(), pd.DataFrame()

        logging.info(f"Found {len(common_patients)} common patients between mRNA and miRNA data.")
        return gene_expr[common_patients], mirna_expr[common_patients]
        
    except FileNotFoundError as e:
        logging.error(f"FATAL: Expression file not found: {e.filename}. Please check paths.")
        return None, None

def load_candidate_interactions(mirtarbase_path, targetscan_path):
    """Tải và hợp nhất các tương tác ứng viên."""
    logging.info("Loading and merging candidate interactions...")
    try:
        validated = pd.read_csv(mirtarbase_path)
        predicted = pd.read_csv(targetscan_path)
        
        validated['validated'] = True
        
        all_interactions = pd.concat([validated, predicted])
        all_interactions = all_interactions.drop_duplicates(subset=['mirna_id', 'gene_id'], keep='first').reset_index(drop=True)
        all_interactions['validated'] = all_interactions['validated'].fillna(False)
        
        logging.info(f"Loaded {len(all_interactions)} unique candidate interactions.")
        return all_interactions
    except FileNotFoundError as e:
        logging.error(f"FATAL: Processed interaction file not found: {e.filename}.")
        return None

def calculate_correlation_weights(gene_expr, mirna_expr, candidates):
    """Tính toán trọng số dựa trên tương quan."""
    if gene_expr.empty or mirna_expr.empty:
        logging.error("Expression data is empty. Cannot calculate correlations.")
        return pd.DataFrame()

    logging.info("Calculating correlation-based weights...")
    
    valid_genes = gene_expr.index
    valid_mirnas = mirna_expr.index
    candidates = candidates[candidates['gene_id'].isin(valid_genes) & candidates['mirna_id'].isin(valid_mirnas)]
    
    results = []
    for _, row in tqdm(candidates.iterrows(), total=candidates.shape[0], desc="Correlating pairs"):
        mirna, gene, is_validated = row['mirna_id'], row['gene_id'], row['validated']
        
        mirna_vec = mirna_expr.loc[mirna].astype(float)
        gene_vec = gene_expr.loc[gene].astype(float)    
        
        if mirna_vec.var() == 0 or gene_vec.var() == 0: continue

        r, p_value = pearsonr(mirna_vec, gene_vec)

        if p_value < P_VALUE_THRESHOLD and r < 0:
            weight = abs(r)
            if is_validated:
                weight = min(1.0, weight + VALIDATED_BONUS)
            
            results.append({'mirna_id': mirna, 'gene_id': gene, 'weight': weight})
            
    return pd.DataFrame(results)

if __name__ == "__main__":
    if not (os.path.exists(MI_RTARBASE_PATH) and os.path.exists(TARGETSCAN_PATH)):
        logging.error("Processed interaction files not found! Please run 'preprocess_interaction_data.py' first.")
    else:
        gene_df, mirna_df = load_expression_data(GENE_EXPR_PATH, MIRNA_EXPR_PATH)
        
        if gene_df is not None and not gene_df.empty:
            candidate_df = load_candidate_interactions(MI_RTARBASE_PATH, TARGETSCAN_PATH)
            
            if candidate_df is not None:
                final_edges_df = calculate_correlation_weights(gene_df, mirna_df, candidate_df)

                if not final_edges_df.empty:
                    output_df = final_edges_df.sort_values(by='weight', ascending=False)
                    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
                    output_df.to_csv(OUTPUT_PATH, index=False)
                    logging.info(f"--- PIPELINE FINISHED ---")
                    logging.info(f"Successfully saved {len(output_df)} edges to {OUTPUT_PATH}")
                else:
                    logging.warning("Pipeline finished, but no significant edges were found with the current thresholds.")