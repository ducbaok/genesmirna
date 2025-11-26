import pandas as pd
import os

# --- CẤU HÌNH ĐƯỜNG DẪN (Chỉnh lại cho đúng máy bạn) ---
GENE_EXPR = 'data/features/genes_expr.txt'  # File biểu hiện Gene
MIRNA_EXPR = 'data/features/mirnas.tsv'     # File biểu hiện miRNA
MIRTARBASE = 'data/processed/mirtarbase_processed.csv'
TARGETSCAN = 'data/processed/targetscan_processed.csv'

def check_overlap():
    print("--- DIAGNOSTIC REPORT ---")
    
    # 1. Kiểm tra miRNA Expression
    try:
        mir_df = pd.read_csv(MIRNA_EXPR, sep='\t', index_col=0)
        print(f"[Expression] miRNAs loaded: {len(mir_df)}")
        print(f"   > Sample IDs (first 3): {list(mir_df.index[:3])}")
    except Exception as e:
        print(f"[Error] Cannot read miRNA file: {e}")
        return

    # 2. Kiểm tra Gene Expression
    try:
        gene_df = pd.read_csv(GENE_EXPR, sep='\t', index_col=0) # Index col có thể là 'Hugo_Symbol'
        print(f"[Expression] Genes loaded: {len(gene_df)}")
        print(f"   > Sample IDs (first 3): {list(gene_df.index[:3])}")
    except Exception as e:
        print(f"[Error] Cannot read Gene file: {e}")
        return

    # 3. Kiểm tra Interaction Data
    for name, path in [('miRTarBase', MIRTARBASE), ('TargetScan', TARGETSCAN)]:
        print(f"\nChecking {name}...")
        try:
            inter_df = pd.read_csv(path)
            unique_mirs = set(inter_df['mirna_id'])
            unique_genes = set(inter_df['gene_id'])
            
            print(f"   > Total interactions: {len(inter_df)}")
            print(f"   > Unique miRNAs in interactions: {len(unique_mirs)}")
            print(f"   > Unique Genes in interactions: {len(unique_genes)}")
            
            # --- KIỂM TRA GIAO NHAU (INTERSECTION) ---
            # Xem bao nhiêu miRNA trong file tương tác CÓ MẶT trong file biểu hiện
            valid_mirs = unique_mirs.intersection(mir_df.index)
            print(f"   > [CRITICAL] miRNAs found in Expression data: {len(valid_mirs)} / {len(unique_mirs)} ({len(valid_mirs)/len(unique_mirs)*100:.1f}%)")
            
            if len(valid_mirs) == 0:
                print("     !!! CẢNH BÁO: Không khớp miRNA ID nào. Kiểm tra lại định dạng (ví dụ: 'hsa-miR-122-5p' vs 'MIMAT...')")
                print(f"     Example Interaction ID: {list(unique_mirs)[0]}")
                print(f"     Example Expression ID:  {list(mir_df.index)[0]}")

            # Xem bao nhiêu Gene trong file tương tác CÓ MẶT trong file biểu hiện
            valid_genes = unique_genes.intersection(gene_df.index)
            print(f"   > [CRITICAL] Genes found in Expression data: {len(valid_genes)} / {len(unique_genes)} ({len(valid_genes)/len(unique_genes)*100:.1f}%)")
            
            if len(valid_genes) == 0:
                print("     !!! CẢNH BÁO: Không khớp Gene ID nào. Kiểm tra lại định dạng (ví dụ: Symbol vs ENSG)")
        except Exception as e:
            print(f"   > Cannot read file {path}: {e}")

if __name__ == "__main__":
    check_overlap()