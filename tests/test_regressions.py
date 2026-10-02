import contextlib
import io
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

os.environ['MPLBACKEND'] = 'Agg'
import matplotlib.pyplot as plt
import pandas as pd
from rmtools import normalize, rm_track, agp_track, depth_track, plot_panel, plot_multi
from rmtools.universal import parse_region


class RegressionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.p = Path(self.tmp.name)
        self.rm = self.p / 'rm.tsv'
        self.df = pd.DataFrame([
            dict(chrom='X', start=100, end=150, strand='+', repeat_name='a',
                 repeat_class='LINE', repeat_family='NA'),
            dict(chrom='X', start=125, end=175, strand='-', repeat_name='b',
                 repeat_class='DNA', repeat_family='hAT')])
        self.df.to_csv(self.rm, sep='\t', index=False)
        self.depth = self.p / 'depth.tsv'
        self.depth.write_text('X\t101\t10\nX\t150\t20\nX\t151\t30\nX\t200\t40\n')
        self.agp = self.p / 'x.agp'
        self.agp.write_text('X\t101\t200\t1\tW\tc1\t1\t100\t+\n')

    def tearDown(self):
        plt.close('all')
        self.tmp.cleanup()

    def cli(self, *args, ok=True):
        proc = subprocess.run([sys.executable, '-m', 'rmtools', *map(str,args)],
                              capture_output=True, text=True)
        if ok:
            self.assertEqual(proc.returncode, 0, proc.stderr)
        else:
            self.assertNotEqual(proc.returncode, 0)
        return proc

    def test_empty_normalize_and_filter(self):
        src = self.p / 'empty.out'; src.write_text('')
        out = self.p / 'empty.tsv'
        self.cli('normalize', '--rm-out', src, '--out', out, '--strain', 's')
        self.assertTrue(rm_track.load_data(out).empty)
        src.write_text('100 0 0 0 X 1 10 (0) C r Simple_repeat 1 10 (0) 1\n')
        self.cli('normalize', '--rm-out', src, '--out', out, '--strain', 's')
        df = rm_track.load_data(out)
        self.assertEqual(df.iloc[0].start, 0)
        self.assertEqual(df.iloc[0].end, 10)
        self.assertEqual(df.iloc[0].strand, '-')
        self.assertEqual(rm_track.choose_taxonomy(df, 'family').iloc[0], 'Simple_repeat/NA')
        self.cli('normalize', '--rm-out', src, '--out', out, '--strain', 's', '--contig', 'missing')
        self.assertTrue(rm_track.load_data(out).empty)

    def test_bins_union_partial_and_empty(self):
        tax = rm_track.choose_taxonomy(self.df, 'class')
        for fun in (rm_track.bin_intervals_dominant, rm_track.bin_intervals_repeat_composition):
            bins = fun(self.df, tax, 50, start=110, end=175)
            totals = bins.groupby(['bin_start','bin_end']).coverage.sum()
            self.assertEqual(totals.to_dict(), {(110,150):40, (150,175):25})
            unann = bins[bins.taxonomy == 'Unannotated']
            self.assertEqual(unann.coverage.sum(), 0)
            empty = fun(self.df.iloc[:0], tax.iloc[:0], 50, start=100, end=200)
            self.assertEqual(empty.coverage.tolist(), [50,50])
        exact = rm_track.bin_intervals_repeat_composition(self.df.iloc[:1], tax.iloc[:1], 50)
        self.assertEqual(exact.bin_start.unique().tolist(), [0,50,100])
        self.assertTrue(rm_track.bin_intervals_repeat_composition(self.df.iloc[:0],tax.iloc[:0],50).empty)

    def test_depth_conversion_and_boundaries(self):
        df = depth_track.load_depth(self.depth)
        self.assertEqual(df.pos.tolist(), [100,149,150,199])
        sub = depth_track.subset_depth(df,'X',100,150)
        self.assertEqual(sub.pos.tolist(),[100,149])
        self.assertEqual(depth_track.bin_depth(sub,50).depth.tolist(),[15])

    def test_agp_width_rebase_and_one_base(self):
        df = agp_track.load_agp(self.agp)
        for rebase,expected in ((True,(0,100)),(False,(100,200))):
            fig,ax=plt.subplots()
            agp_track.plot_agp_layers(df,ax,region_start=100,region_end=200,rebase=rebase)
            x=ax.collections[0].get_paths()[0].vertices[:,0]
            self.assertEqual((x.min(),x.max()),expected)
        self.agp.write_text('X\t1\t1\t1\tW\tc1\t1\t1\t+\n')
        df=agp_track.load_agp(self.agp)
        self.assertEqual(len(agp_track.subset_agp(df,'X',0,1)),1)
        fig,ax=plt.subplots();agp_track.plot_agp_layers(df,ax)
        x=ax.collections[0].get_paths()[0].vertices[:,0]
        self.assertEqual((x.min(),x.max()),(0,1))

    def test_panel_alignment_and_empty(self):
        axes=plot_panel.plot_panel('X:100-200',rm_path=self.rm,depth_path=self.depth,
                                  agp_path=self.agp,rm_bin_size=50,depth_bin_size=50)
        self.assertTrue(all(ax.get_xlim()==(100,200) for ax in axes))
        self.assertEqual(axes[1].lines[0].get_xdata().tolist(),[100,150])
        self.assertTrue(all(b.get_x()>=100 for b in axes[0].patches))
        axes=plot_panel.plot_panel('X:300-400',rm_path=self.rm,rm_bin_size=50)
        self.assertEqual(sum(b.get_height() for b in axes[0].patches),100)

    def test_multi_colors_origin_and_empty(self):
        second=self.p/'second.tsv'
        self.df.iloc[:1].to_csv(second,sep='\t',index=False)
        control=pd.DataFrame({'path':[str(self.rm),str(second)],'contig':['X','X'],'label':['a','b']})
        fig=plot_multi.plot_multi(control,'class',50)
        cmap=rm_track.make_color_map(['DNA','LINE'])
        for ax in fig.axes:
            self.assertTrue(any(b.get_facecolor()==cmap['LINE'] for b in ax.patches))
            self.assertTrue(any(b.get_x()==0 and b.get_height()==50 for b in ax.patches))
        control['contig']='missing'
        fig=plot_multi.plot_multi(control,'class',50)
        self.assertEqual(len(fig.axes),2)

    def test_plot_commands(self):
        sizes=self.p/'sizes.tsv';sizes.write_text('X\t300\nmissing\t300\n')
        for region in ('X:100-200','missing:0-100'):
            for cmd,flag,path in (('plot-contig','--rm',self.rm),('depth-track','--depth',self.depth),('agp-track','--agp',self.agp)):
                self.cli(cmd,flag,path,'--region',region,'--out',self.p/(cmd+'.pdf'))
        self.cli('plot-contig','--rm',self.rm,'--region','X','--taxonomy','family','--bin-size','50','--sizes',sizes,'--out',self.p/'contig.pdf')
        for cmd in ('plot-main','plot-assembly'):
            self.cli(cmd,'--rm',self.rm,'--main','X','missing','--sizes',sizes,'--out',self.p/(cmd+'.pdf'))
        self.cli('panel','--rm',self.rm,'--depth',self.depth,'--agp',self.agp,'--region','X:100-200','--out',self.p/'panel.pdf')
        control=self.p/'multi.tsv';control.write_text(f'path\tcontig\tlabel\ns{0}\tX\ta\n'.replace('s0',str(self.rm)))
        self.cli('plot-multi','--control',control,'--bin-size','50','--out',self.p/'multi.pdf')
        control=self.p/'compare.tsv'
        control.write_text(f'strain\trm\tdepth\tagp\ns1\t{self.rm}\t\t\ns2\t\t{self.depth}\t{self.agp}\n')
        self.cli('plot-compare','--control',control,'--contigs','X','missing','--out-prefix',self.p/'compare')

    def test_invalid_inputs(self):
        for region in ('X:20-10','X:1-1','X:-1-10',':0-10'):
            with self.assertRaises(ValueError):parse_region(region)
        self.cli('plot-contig','--rm',self.rm,'--region','X','--bin-size','0','--out',self.p/'bad.pdf',ok=False)
        with self.assertRaises(ValueError):
            rm_track.bin_intervals_dominant(self.df,rm_track.choose_taxonomy(self.df,'class'),0)


if __name__ == '__main__':
    unittest.main()
