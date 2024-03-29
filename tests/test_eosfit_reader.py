# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c), Gianfranco Ulian and Giovanni Valdre'.                      #
# All rights reserved.                                                       #
#                                                                            #
# This file is part of the Quantas code.                                     #
#                                                                            #
# For further information on the license, see the LICENSE file               #
##############################################################################

""" Test cases for the input file readers. """

import os
import unittest
import tempfile
import numpy as np

from quantas.IO.eos_reader import EoSInputFileReader


class EoSReaderTest(unittest.TestCase):

    valid_eos_input = '''JOB
Topaz (Gatta et al., 2014)
FORMAT
P	sigmaP	V	sigmaV	a	sigmaA	b	sigmaB	c	sigmaC
DATA
0.0001	1.e-7	345.46	0.02	4.6627	0.0002	8.8343	0.0004	8.3867	0.0002
0.41	0.05	344.80	0.02	4.6597	0.0002	8.8303	0.0003	8.3799	0.0002
1.88	0.05	341.78	0.02	4.6451	0.0002	8.8108	0.0003	8.3509	0.0002
2.57	0.05	340.35	0.02	4.6378	0.0002	8.8022	0.0003	8.3372	0.0002
3.14	0.05	339.15	0.02	4.6325	0.0002	8.7938	0.0003	8.3254	0.0002
3.67	0.05	338.05	0.02	4.6267	0.0002	8.7872	0.0003	8.3150	0.0002
4.65	0.05	336.28	0.02	4.6179	0.0002	8.7759	0.0003	8.2976	0.0002
5.16	0.05	335.34	0.02	4.6135	0.0002	8.7692	0.0003	8.2887	0.0002
5.79	0.05	334.08	0.02	4.6071	0.0002	8.7614	0.0003	8.2765	0.0002
6.52	0.05	332.82	0.02	4.6009	0.0002	8.7533	0.0003	8.2641	0.0002
6.62	0.10	331.40	0.20	4.579	0.001	8.751	0.002	8.269	0.001
7.41	0.05	331.23	0.03	4.5930	0.0004	8.7427	0.0005	8.2488	0.0003
8.37	0.05	329.45	0.03	4.5841	0.0004	8.7308	0.0005	8.2316	0.0003
8.97	0.05	328.53	0.04	4.5795	0.0005	8.7247	0.0006	8.2224	0.0004
9.38	0.10	326.80	0.20	4.568	0.001	8.714	0.002	8.209	0.002
11.34	0.10	324.10	0.20	4.559	0.001	8.690	0.002	8.180	0.002
15.27	0.10	316.10	0.20	4.525	0.001	8.636	0.002	8.089	0.001
19.81	0.10	309.90	0.20	4.493	0.001	8.589	0.002	8.029	0.002
26.02	0.10	300.00	0.20	4.446	0.001	8.512	0.002	7.927	0.002
31.49	0.10	294.00	0.20	4.414	0.001	8.464	0.002	7.869	0.002
37.82	0.10	289.20	0.20	4.399	0.001	8.412	0.002	7.814	0.002
41.00	0.10	284.00	0.20	4.364	0.001	8.377	0.002	7.767	0.002
43.84	0.10	283.00	0.20	4.363	0.001	8.356	0.002	7.763	0.002
45.26	0.10	281.60	0.20	4.352	0.001	8.346	0.002	7.753	0.002
45.33	0.10	281.70	0.20	4.360	0.001	8.348	0.003	7.739	0.003
'''
        
    def test_reading(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
            filename = tmp.name
            basename = os.path.split(filename)[1]
            tmp.write(self.valid_eos_input)
            tmp.flush()
            reader = EoSInputFileReader(filename)

        self.assertTrue(reader.completed, 'Reading èrocess should be completed')
        self.assertEqual(reader.jobname, 'Topaz (Gatta et al., 2014)',
                         "Jobname should be 'Topaz (Gatta et al., 2014)'")
        self.assertTrue(isinstance(reader.data['p'], np.ndarray),
                        'Pressure data should be present')
        self.assertTrue(isinstance(reader.data['sigmap'], np.ndarray),
                        'Uncertainties on pressure data should be present')
        
        self.assertTrue(isinstance(reader.data['v'], np.ndarray),
                        'Volume data should be present')
        self.assertTrue(isinstance(reader.data['sigmav'], np.ndarray),
                        'Uncertainties on volume data should be present')
        
        self.assertTrue(isinstance(reader.data['a'], np.ndarray),
                        'a-axis data should be present')
        self.assertTrue(isinstance(reader.data['sigmaa'], np.ndarray),
                        'Uncertainties on a-axis data should be present')
        
        self.assertTrue(isinstance(reader.data['b'], np.ndarray),
                        'b-axis data should be present')
        self.assertTrue(isinstance(reader.data['sigmab'], np.ndarray),
                        'Uncertainties on b-axis data should be present')
        
        self.assertTrue(isinstance(reader.data['c'], np.ndarray),
                        'c-axis data should be present')
        self.assertTrue(isinstance(reader.data['sigmac'], np.ndarray),
                        'Uncertainties on c-axis data should be present')
            
        tmp.close()
        os.unlink(tmp.name)


if __name__ == '__main__':
    unittest.main()
