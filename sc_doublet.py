import os
import glob
import scanpy as sc
import scrublet as scr
import scipy.sparse as sp

def get_sample_id(h5_file_path):
    filename = os.path.basename(h5_file_path)
    sample_id = "_".join(filename.split("_")[0:3])  # Extract CG_SB_NBxxxxxxx
    return sample_id

def run_scrublet(adata_x):
    try:
        counts_matrix = adata_x.X
        if not sp.issparse(counts_matrix):
            counts_matrix = sp.csr_matrix(counts_matrix)

        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.1, sim_doublet_ratio=3)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()

        # Now apply a manual threshold on the scores:
        threshold = 0.25
        predicted_doublets_strict = doublet_scores > threshold

        adata_x.obs["scrublet"] = ["doublet" if x else "singlet" for x in predicted_doublets_strict]
        adata_x.obs["scrublet_score"] = doublet_scores

    except Exception as e:
        print(f"Scrublet failed: {e}")

    return adata_x

def load_and_merge_data(h5_file_list):
    adatas = []
    for h5_file in h5_file_list:
        sample_id = get_sample_id(h5_file)

        adata = sc.read_10x_h5(h5_file)
        adata.var_names_make_unique()

        adata = run_scrublet(adata)

        adata.obs["CellID"] = [sample_id + ":" + bc for bc in adata.obs.index]
        adata.obs["SampleID"] = sample_id
        adata.obs["Barcode"] = adata.obs.index

        adatas.append(adata)

    if len(adatas) == 1:
        return adatas[0]
    else:
        adata_merged = adatas[0].concatenate(*adatas[1:], index_unique="-")
        return adata_merged

def main():
    base_dir = "/nfs/users/nfs_l/lr26/"
    h5_files = glob.glob(os.path.join(base_dir, "*filtered_feature_bc_matrix.h5"))
    
    print(f"Found {len(h5_files)} files")

    adata_merged = load_and_merge_data(h5_files)
    print("Finished loading, Scrublet, and merging data")

    out_path = os.path.join(base_dir, "adata_merged_with_scrublet.h5ad")
    adata_merged.write(out_path)
    print(f"Saved merged AnnData with Scrublet scores to {out_path}")

if __name__ == "__main__":
    main()
