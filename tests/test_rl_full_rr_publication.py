"""Read-only presentation contracts; synthetic tables, no spectral solves."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.analysis import plot_rl_full_rr_publication as figures


class PublicationTables(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.spec = dict(name="example", mode="BETA", fixed=dict(mu=.1, tau=0), grid=[0., .5, .6, 1.])
        self.config = dict(maps=[self.spec], material=dict(epsilon=.005),
            grid_provenance=dict(analysis_step=.5, local_step=.1, local_windows=dict(example=[.5, 1.])))
        self.rows = [dict(map="example", x=x, rank=k, status="PASS", Lambda=k+x/10)
                     for x in self.spec["grid"] for k in range(1,10)]
        self.write()

    def write(self):
        table = self.folder/"spectrum.csv"
        with table.open("w", newline="") as f:
            w=csv.DictWriter(f,fieldnames=self.rows[0]); w.writeheader(); w.writerows(self.rows)
        manifest=dict(status="PASS", artifacts={"spectrum.csv":figures.archive.sha(table)},
            identity=dict(maps=[self.spec],material=self.config["material"],
                          local_policies=dict(example=dict(spectrum_semantics="sorted_positions"))))
        (self.folder/"manifest.json").write_text(json.dumps(manifest))

    def test_altered_table_and_wrong_configuration_rejected(self):
        with (self.folder/"spectrum.csv").open("a") as f: f.write("\n")
        with self.assertRaisesRegex(ValueError,"checksum"):
            figures.read_maps(self.folder,self.config)
        self.write()
        wrong=dict(self.config,material=dict(epsilon=.01))
        with self.assertRaisesRegex(ValueError,"Configuration"):
            figures.read_maps(self.folder,wrong)

    def test_missing_guard_duplicate_and_failure_rejected(self):
        self.rows.pop(); self.write()
        with self.assertRaises(KeyError): figures.read_maps(self.folder,self.config)
        self.rows.append(self.rows[0]); self.write()
        with self.assertRaisesRegex(ValueError,"Duplicate"): figures.read_maps(self.folder,self.config)
        self.rows.pop(); self.rows[0]["status"]="UNRESOLVED"; self.write()
        with self.assertRaisesRegex(ValueError,"Unresolved"): figures.read_maps(self.folder,self.config)

    def test_curves_keep_saved_nodes_and_hide_guard(self):
        _,maps=figures.read_maps(self.folder,self.config)
        fig,ax=plt.subplots(); self.addCleanup(plt.close,fig)
        figures.draw(ax,maps[0],range(1,9),ymax=10)
        self.assertEqual(len(ax.lines),8)
        for k,line in enumerate(ax.lines,1):
            self.assertEqual(list(line.get_xdata()),[0.,.5,.6,1.])
            self.assertEqual(list(line.get_ydata()),[k+x/10 for x in [0.,.5,.6,1.]])
        maps[0]['y'][2][1]=maps[0]['y'][2][0]+.01
        gap=figures.sampled_gaps(maps,self.config)[0]
        self.assertEqual(gap['beta_deg'],.6)
        self.assertEqual(gap['grid'],'local')


if __name__ == "__main__":
    unittest.main()
