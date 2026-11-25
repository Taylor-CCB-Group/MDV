import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

from mdvtools.conversions import convert_mudata_to_mdv


def _mudata_with_pulled_obs():
    n = 4
    obs_names = pd.Index([f"c{i}" for i in range(n)])
    batch = ["b1", "b1", "b2", "b2"]
    rna = ad.AnnData(
        np.ones((n, 3), dtype=np.float32),
        obs=pd.DataFrame({"condition": ["control"] * n, "batch": batch}, index=obs_names),
    )
    prot = ad.AnnData(
        np.ones((n, 2), dtype=np.float32),
        obs=pd.DataFrame({"batch": batch}, index=obs_names),
    )
    md = mu.MuData({"rna": rna, "prot": prot})
    md.pull_obs()  # adds rna:condition, rna:batch, prot:batch
    md.obs["condition"] = ["treated"] * n
    md.obs["batch"] = batch
    return md


def test_mudata_drops_only_identical_modality_copies(tmp_path):
    project = convert_mudata_to_mdv(str(tmp_path / "proj"), _mudata_with_pulled_obs(), delete_existing=True)

    fields = {c["field"] for c in project.get_datasource_metadata("cells")["columns"]}
    # same values as the global column: redundant
    assert "rna:batch" not in fields
    assert "prot:batch" not in fields
    # differs from the global column: a separate annotation that must survive
    assert "rna:condition" in fields
    assert list(project.get_column("cells", "rna:condition")) == ["control"] * 4
    assert list(project.get_column("cells", "condition")) == ["treated"] * 4
