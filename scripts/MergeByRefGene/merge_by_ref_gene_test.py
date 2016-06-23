from unittest.case import TestCase
from StringIO import StringIO
from merge_by_ref_gene import merge_by_ref_gene


class MergeByRefGeneTest(TestCase):
    def setUp(self):
        self.addTypeEqualityFunc(str, self.assertMultiLineEqual)

    def test_simple(self):
        aminos = StringIO("""\
rname,gene,aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
r1,NS3,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
r1,NS3,1,0,1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0
r1,NS3,2,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0
r1,NS3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0
r2,NS3,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
r2,NS3,1,0,0,2,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0
r2,NS3,2,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0
r2,NS3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0
""")
        expected_merged = """\
gene,ref.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
NS3,1,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
NS3,2,0,1,2,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0
NS3,3,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0
NS3,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0
"""
        merged = StringIO()

        merge_by_ref_gene(aminos, merged)

        self.maxDiff = None
        self.assertEqual(expected_merged, merged.getvalue())

    def test_consensus(self):
        """ Make sure 10 A's are chosen over 2 S's to make APIT, not SPIT. """
        aminos = StringIO("""\
rname,gene,aa.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
r1,NS3,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0
r1,NS3,1,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0,0,0,0,0
r1,NS3,2,0,0,0,0,0,0,0,12,0,0,0,0,0,0,0,0,0,0,0,0,0
r1,NS3,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0
""")
        expected_merged = """\
gene,ref.pos,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,*
NS3,1,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0
NS3,2,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0,0,0,0,0
NS3,3,0,0,0,0,0,0,0,12,0,0,0,0,0,0,0,0,0,0,0,0,0
NS3,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,0,0,0,0
"""
        merged = StringIO()

        merge_by_ref_gene(aminos, merged)

        self.maxDiff = None
        self.assertEqual(expected_merged, merged.getvalue())
