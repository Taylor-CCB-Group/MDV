from pathlib import Path
import h5py
import numpy as np
import scipy.sparse
import scanpy as sc
from scanpy import AnnData
from typing import cast
import json
import traceback

def run(workspace: str) -> None:
    ws = Path(workspace)
    try:
        with h5py.File(ws / "input" / "tray.h5", "r") as f:
            x = cast(h5py.Dataset, f["x"])[:]
            i = cast(h5py.Dataset, f["i"])[:]
            p = cast(h5py.Dataset, f["p"])[:]
            n_cells = int(cast(int, f.attrs["n_cells"]))
            n_genes = int(cast(int, f.attrs["n_genes"]))
            out_name = cast(str, f.attrs["output_name"])
            neighbors_kwargs = json.loads(cast(str, f.attrs.get("kwargs.neighbors", "{}")))
            umap_kwargs = json.loads(cast(str, f.attrs.get("kwargs.umap", "{}")))
            pca_kwargs = json.loads(cast(str, f.attrs.get("kwargs.pca", "{}")))

        # rebuild the CSC matrix from the tray triplet - n_cells from metadata, NOT max(i) + 1
        X = scipy.sparse.csc_matrix((x, i, p), shape=(n_cells, n_genes))

        # inline sanitise (inspired from _coerce_x_to_csc_without_nonfinite)
        X.data = np.nan_to_num(X.data, nan=0.0, posinf=0.0, neginf=0.0)
        X.eliminate_zeros()

        # PCA -> neighbour graph on X_pca -> embed; a graph on every gene lets noise swamp the clusters
        adata = AnnData(X=X)
        if "n_comps" in pca_kwargs:
            # arpack needs fewer components than the smaller dimension; manifest records what was used
            pca_kwargs["n_comps"] = min(int(pca_kwargs["n_comps"]), min(n_cells, n_genes) - 1)
        sc.pp.pca(adata, **pca_kwargs)
        n_comps = int(adata.obsm["X_pca"].shape[1])   # components actually used
        sc.pp.neighbors(adata, use_rep="X_pca", **neighbors_kwargs)
        sc.tl.umap(adata, **umap_kwargs)
        embedding = np.asarray(adata.obsm["X_umap"], dtype=np.float64) # (n_cells, n_components)

        with h5py.File(ws / "output" / "result.h5", "w") as f:
            for dim in range(embedding.shape[1]):
                # 1-indexed {base}_{i+1}, mirroring _add_dims: UMAP -> UMAP_1, UMAP_2
                f.create_dataset(f"{out_name}_{dim + 1}", data=embedding[:, dim])
            f.attrs["output_name"] = out_name

        (ws / "output" / "manifest.json").write_text(
            json.dumps({"rows": int(embedding.shape[0]), "columns": int(embedding.shape[1]), "n_comps": n_comps})
        ) # basic provenance
        (ws / "STATUS").write_text("done")

    except Exception:
        (ws / "STATUS").write_text("failed")
        (ws / "output" / "error.txt").write_text(traceback.format_exc(), encoding="utf-8", errors="replace")
        raise
