"""Figure checks only: no new eigenproblem, tracking or recovery."""
import sys
import numpy as np
import pytest
from scripts.lib import inplane_spring_figure03 as plot


@pytest.fixture(scope='module')
def saved():
    if not (plot.SOURCE/'shapes.npz').exists():
        pytest.skip('saved K11 shapes required; do not reconstruct automatically')
    return plot.load_saved()


def test_fixed_pair_and_no_missing_connections(saved):
    assert len(saved['data'])==44 and len(saved['angles'])==11
    assert saved['snapshots']==[25.340576171875,30.340576171875,35.340576171875]
    assert saved['event']['symmetric_candidate']['classification']=='CROSSING_SUPPORTED'
    assert saved['event']['classification']=='AVOIDED_CROSSING_RESOLVED'
    for mu in (0.,.01):
        for branch in plot.BRANCHES:
            for b in saved['angles']:
                r=saved['data'][(mu,branch,b)]['row']
                assert r['root_status']==r['mapping_status']=='CONFIRMED'
                assert r['branch_id']==branch and r['model']=='EB'


def test_shapes_change_only_by_one_global_sign(saved):
    with np.load(plot.SOURCE/'shapes.npz',allow_pickle=False) as source:
        for p in saved['phase_records']:
            item=saved['data'][(p['mu'],p['branch_id'],p['beta_deg'])]
            np.testing.assert_array_equal(item['states'],p['sign']*source[item['row']['shape_key']+'__states'])
    assert plot.COLORS=={'mode_01':'#0072B2','mode_02':'#D55E00'}


@pytest.mark.parametrize('mu',[0.,.01])
def test_physical_coordinates_and_common_display_scale(saved,mu):
    for branch in plot.BRANCHES:
        for beta in saved['snapshots']:
            item=saved['data'][(mu,branch,beta)]
            reference,disp=plot.global_centrelines(beta,mu,item['states'])
            np.testing.assert_allclose(np.linalg.norm(reference[:,0],axis=1),[1-mu,1+mu],rtol=1e-14)
            np.testing.assert_array_equal(reference[:,-1],np.zeros((2,2)))
            np.testing.assert_allclose(disp[:,0],np.zeros((2,2)),atol=1e-13)
            np.testing.assert_allclose(disp[0,-1],disp[1,-1],atol=2e-9,rtol=0.)
            # Orthogonal display rotation preserves displacement norm and arm amplitudes.
            np.testing.assert_allclose(np.linalg.norm(disp,axis=2),np.linalg.norm(item['states'][:,:,:2],axis=2),atol=1e-13)
            assert np.linalg.norm(saved['display_scale']*disp,axis=2).max()<=plot.DISPLAY_AMPLITUDE+1e-14


def test_plotted_vertices_are_saved_roots(saved):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for mu in (0.,.01):
        fig=plot.draw_case(saved,mu); ax=fig.axes[0]
        lines={line.get_label():line for line in ax.lines if line.get_label() in ('1','2')}
        for i,branch in enumerate(plot.BRANCHES,1):
            np.testing.assert_array_equal(lines[str(i)].get_xdata(),saved['angles'])
            np.testing.assert_array_equal(lines[str(i)].get_ydata(),
                [float(saved['data'][(mu,branch,b)]['row']['Lambda']) for b in saved['angles']])
        assert len(fig.axes)==7
        plt.close(fig)


def test_existing_entry_plot_only_no_science(saved,tmp_path,monkeypatch):
    from scripts.analysis.laminated_beams import check_inplane_spring_robustness as run
    def forbidden(*a,**k):pytest.fail('figure preset must not compute roots, matrices, forms or tracking')
    for obj,name in ((run,'Run'),(run,'preflight'),(run,'track_window'),(run.mechanics,'recover'),
                     (run.eb,'state_matrix'),(run.eb,'transfer_matrix'),
                     (run.rlb,'state_matrix'),(run.rlb,'transfer_matrix')):
        monkeypatch.setattr(obj,name,forbidden)
    before={p:plot.sha(plot.SOURCE/p) for p in plot.SOURCE_FILES}
    original=plot.render
    monkeypatch.setattr(plot,'render',lambda:original(output=tmp_path))
    monkeypatch.setattr(sys,'argv',['check_inplane_spring_robustness.py','plot-only','--figure03'])
    run.main()
    assert before=={p:plot.sha(plot.SOURCE/p) for p in plot.SOURCE_FILES}
    manifest=plot.read_json(tmp_path/'figure_manifest.json')
    assert all(manifest[k]==0 for k in ('new_roots','new_beta','matrix_calls','shape_recoveries','tracking_calls','interpolated_frequencies'))
    from PIL import Image
    for stem in ('figure03a_crossing','figure03b_veering'):
        with Image.open(tmp_path/(stem+'.png')) as im:
            assert im.size==(2820,2160) and abs(im.info['dpi'][0]-300)<.1
        pdf=(tmp_path/(stem+'.pdf')).read_bytes()
        assert pdf.startswith(b'%PDF-') and b'/Subtype /Image' not in pdf  # fully vector figure
