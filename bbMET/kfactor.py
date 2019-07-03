import os, sys

#//--------------------------------------------------------------------------------------
def getEWKW(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 0.94269
  elif(pt>=170  and pt < 200)  : weight = 0.902615
  elif(pt>=200  and pt < 230)  : weight = 0.898827
  elif(pt>=230  and pt < 260)  : weight = 0.959081
  elif(pt>=260  and pt < 290)  : weight = 0.891248
  elif(pt>=290  and pt < 320)  : weight = 0.860188
  elif(pt>=320  and pt < 350)  : weight = 0.884811
  elif(pt>=350  and pt < 390)  : weight = 0.868131
  elif(pt>=390  and pt < 430)  : weight = 0.848655
  elif(pt>=430  and pt < 470)  : weight = 0.806186
  elif(pt>=470  and pt < 510)  : weight = 0.848507
  elif(pt>=510  and pt < 550)  : weight = 0.83763
  elif(pt>=550  and pt < 590)  : weight = 0.792152
  elif(pt>=590  and pt < 640)  : weight = 0.730731
  elif(pt>=640  and pt < 690)  : weight = 0.778061
  elif(pt>=690  and pt < 740)  : weight = 0.771811
  elif(pt>=740  and pt < 790)  : weight = 0.795004
  elif(pt>=790  and pt < 840)  : weight = 0.757859
  elif(pt>=840  and pt < 900)  : weight = 0.709571
  elif(pt>=900  and pt < 960)  : weight = 0.702751
  elif(pt>=960  and pt < 1020) : weight = 0.657821
  elif(pt>=1020 and pt < 1090) : weight = 0.762559
  elif(pt>=1090 and pt < 1160) : weight = 0.845925
  elif(pt>=1160)              : weight = 0.674034

  return weight


def getQCDW(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 1.43896
  elif(pt>=170  and pt < 200)  : weight = 1.45307
  elif(pt>=200  and pt < 230)  : weight = 1.41551
  elif(pt>=230  and pt < 260)  : weight = 1.42199
  elif(pt>=260  and pt < 290)  : weight = 1.3477
  elif(pt>=290  and pt < 320)  : weight = 1.35302
  elif(pt>=320  and pt < 350)  : weight = 1.34289
  elif(pt>=350  and pt < 390)  : weight = 1.32474
  elif(pt>=390  and pt < 430)  : weight = 1.23267
  elif(pt>=430  and pt < 470)  : weight = 1.22641
  elif(pt>=470  and pt < 510)  : weight = 1.23149
  elif(pt>=510  and pt < 550)  : weight = 1.21593
  elif(pt>=550  and pt < 590)  : weight = 1.16506
  elif(pt>=590  and pt < 640)  : weight = 1.01718
  elif(pt>=640  and pt < 690)  : weight = 1.01575
  elif(pt>=690  and pt < 740)  : weight = 1.05425
  elif(pt>=740  and pt < 790)  : weight = 1.05992
  elif(pt>=790  and pt < 840)  : weight = 1.01503
  elif(pt>=840  and pt < 900)  : weight = 1.01761
  elif(pt>=900  and pt < 960)  : weight = 0.947194
  elif(pt>=960  and pt < 1020) : weight = 0.932754
  elif(pt>=1020 and pt < 1090) : weight = 1.00849
  elif(pt>=1090 and pt < 1160) : weight = 0.94805
  elif(pt>=1160)              : weight = 0.86956

  return weight


def getRenUpW(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 1.33939
  elif(pt>=170  and pt < 200)  : weight = 1.34649
  elif(pt>=200  and pt < 230)  : weight = 1.3115
  elif(pt>=230  and pt < 260)  : weight = 1.30782
  elif(pt>=260  and pt < 290)  : weight = 1.24619
  elif(pt>=290  and pt < 320)  : weight = 1.24531
  elif(pt>=320  and pt < 350)  : weight = 1.23472
  elif(pt>=350  and pt < 390)  : weight = 1.20709
  elif(pt>=390  and pt < 430)  : weight = 1.13368
  elif(pt>=430  and pt < 470)  : weight = 1.12227
  elif(pt>=470  and pt < 510)  : weight = 1.12612
  elif(pt>=510  and pt < 550)  : weight = 1.1118
  elif(pt>=550  and pt < 590)  : weight = 1.06549
  elif(pt>=590  and pt < 640)  : weight = 0.931838
  elif(pt>=640  and pt < 690)  : weight = 0.929282
  elif(pt>=690  and pt < 740)  : weight = 0.959553
  elif(pt>=740  and pt < 790)  : weight = 0.955823
  elif(pt>=790  and pt < 840)  : weight = 0.920614
  elif(pt>=840  and pt < 900)  : weight = 0.917243
  elif(pt>=900  and pt < 960)  : weight = 0.855649
  elif(pt>=960  and pt < 1020) : weight = 0.84587
  elif(pt>=1020 and pt < 1090) : weight = 0.906862
  elif(pt>=1090 and pt < 1160) : weight = 0.858763
  elif(pt>=1160)              : weight = 0.794909

  return weight


def getRenDownW(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 1.52924
  elif(pt>=170  and pt < 200)  : weight = 1.55189
  elif(pt>=200  and pt < 230)  : weight = 1.51041
  elif(pt>=230  and pt < 260)  : weight = 1.53402
  elif(pt>=260  and pt < 290)  : weight = 1.44328
  elif(pt>=290  and pt < 320)  : weight = 1.45859
  elif(pt>=320  and pt < 350)  : weight = 1.4486
  elif(pt>=350  and pt < 390)  : weight = 1.44818
  elif(pt>=390  and pt < 430)  : weight = 1.33058
  elif(pt>=430  and pt < 470)  : weight = 1.33239
  elif(pt>=470  and pt < 510)  : weight = 1.33976
  elif(pt>=510  and pt < 550)  : weight = 1.32336
  elif(pt>=550  and pt < 590)  : weight = 1.26854
  elif(pt>=590  and pt < 640)  : weight = 1.10489
  elif(pt>=640  and pt < 690)  : weight = 1.1055
  elif(pt>=690  and pt < 740)  : weight = 1.15559
  elif(pt>=740  and pt < 790)  : weight = 1.17682
  elif(pt>=790  and pt < 840)  : weight = 1.118
  elif(pt>=840  and pt < 900)  : weight = 1.13097
  elif(pt>=900  and pt < 960)  : weight = 1.04988
  elif(pt>=960  and pt < 1020) : weight = 1.02796
  elif(pt>=1020 and pt < 1090) : weight = 1.12438
  elif(pt>=1090 and pt < 1160) : weight = 1.04704
  elif(pt>=1160)              : weight = 0.94791

  return weight

def getFacUpW(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 1.44488
  elif(pt>=170  and pt < 200)  : weight = 1.45847
  elif(pt>=200  and pt < 230)  : weight = 1.41004
  elif(pt>=230  and pt < 260)  : weight = 1.39864
  elif(pt>=260  and pt < 290)  : weight = 1.3391
  elif(pt>=290  and pt < 320)  : weight = 1.33663
  elif(pt>=320  and pt < 350)  : weight = 1.32073
  elif(pt>=350  and pt < 390)  : weight = 1.3076
  elif(pt>=390  and pt < 430)  : weight = 1.20904
  elif(pt>=430  and pt < 470)  : weight = 1.20066
  elif(pt>=470  and pt < 510)  : weight = 1.20462
  elif(pt>=510  and pt < 550)  : weight = 1.18286
  elif(pt>=550  and pt < 590)  : weight = 1.12586
  elif(pt>=590  and pt < 640)  : weight = 0.990615
  elif(pt>=640  and pt < 690)  : weight = 0.984473
  elif(pt>=690  and pt < 740)  : weight = 1.0171
  elif(pt>=740  and pt < 790)  : weight = 1.01706
  elif(pt>=790  and pt < 840)  : weight = 0.986107
  elif(pt>=840  and pt < 900)  : weight = 0.972452
  elif(pt>=900  and pt < 960)  : weight = 0.910183
  elif(pt>=960  and pt < 1020) : weight = 0.885284
  elif(pt>=1020 and pt < 1090) : weight = 0.950662
  elif(pt>=1090 and pt < 1160) : weight = 0.89605
  elif(pt>=1160)              : weight = 0.823212

  return weight


def getFacDownW(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 1.43694
  elif(pt>=170  and pt < 200)  : weight = 1.45181
  elif(pt>=200  and pt < 230)  : weight = 1.42649
  elif(pt>=230  and pt < 260)  : weight = 1.4513
  elif(pt>=260  and pt < 290)  : weight = 1.3624
  elif(pt>=290  and pt < 320)  : weight = 1.37472
  elif(pt>=320  and pt < 350)  : weight = 1.3718
  elif(pt>=350  and pt < 390)  : weight = 1.34637
  elif(pt>=390  and pt < 430)  : weight = 1.26276
  elif(pt>=430  and pt < 470)  : weight = 1.25762
  elif(pt>=470  and pt < 510)  : weight = 1.26389
  elif(pt>=510  and pt < 550)  : weight = 1.25504
  elif(pt>=550  and pt < 590)  : weight = 1.2116
  elif(pt>=590  and pt < 640)  : weight = 1.04846
  elif(pt>=640  and pt < 690)  : weight = 1.052
  elif(pt>=690  and pt < 740)  : weight = 1.09788
  elif(pt>=740  and pt < 790)  : weight = 1.10967
  elif(pt>=790  and pt < 840)  : weight = 1.046
  elif(pt>=840  and pt < 900)  : weight = 1.06961
  elif(pt>=900  and pt < 960)  : weight = 0.989001
  elif(pt>=960  and pt < 1020) : weight = 0.987722
  elif(pt>=1020 and pt < 1090) : weight = 1.07552
  elif(pt>=1090 and pt < 1160) : weight = 1.00802
  elif(pt>=1160)              : weight = 0.922172

  return weight


#//--------------------------------------------------------------------------------------
def getEWKZ(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 0.970592
  elif(pt>=170  and pt < 200)  : weight = 0.964424
  elif(pt>=200  and pt < 230)  : weight = 0.956695
  elif(pt>=230  and pt < 260)  : weight = 0.948747
  elif(pt>=260  and pt < 290)  : weight = 0.941761
  elif(pt>=290  and pt < 320)  : weight = 0.934246
  elif(pt>=320  and pt < 350)  : weight = 0.927089
  elif(pt>=350  and pt < 390)  : weight = 0.919181
  elif(pt>=390  and pt < 430)  : weight = 0.909926
  elif(pt>=430  and pt < 470)  : weight = 0.900911
  elif(pt>=470  and pt < 510)  : weight = 0.892561
  elif(pt>=510  and pt < 550)  : weight = 0.884353
  elif(pt>=550  and pt < 590)  : weight = 0.8761
  elif(pt>=590  and pt < 640)  : weight = 0.867687
  elif(pt>=640  and pt < 690)  : weight = 0.858047
  elif(pt>=690  and pt < 740)  : weight = 0.849014
  elif(pt>=740  and pt < 790)  : weight = 0.840317
  elif(pt>=790  and pt < 840)  : weight = 0.832017
  elif(pt>=840  and pt < 900)  : weight = 0.823545
  elif(pt>=900  and pt < 960)  : weight = 0.814596
  elif(pt>=960  and pt < 1020) : weight = 0.806229
  elif(pt>=1020 and pt < 1090) : weight = 0.798038
  elif(pt>=1090 and pt < 1160) : weight = 0.789694
  elif(pt>=1160)              : weight = 0.781163

  return weight


def getQCDZ(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 1.47528
  elif(pt>=170  and pt < 200)  : weight = 1.5428
  elif(pt>=200  and pt < 230)  : weight = 1.49376
  elif(pt>=230  and pt < 260)  : weight = 1.39119
  elif(pt>=260  and pt < 290)  : weight = 1.40538
  elif(pt>=290  and pt < 320)  : weight = 1.44661
  elif(pt>=320  and pt < 350)  : weight = 1.38176
  elif(pt>=350  and pt < 390)  : weight = 1.37381
  elif(pt>=390  and pt < 430)  : weight = 1.29145
  elif(pt>=430  and pt < 470)  : weight = 1.33452
  elif(pt>=470  and pt < 510)  : weight = 1.25765
  elif(pt>=510  and pt < 550)  : weight = 1.24265
  elif(pt>=550  and pt < 590)  : weight = 1.24331
  elif(pt>=590  and pt < 640)  : weight = 1.16187
  elif(pt>=640  and pt < 690)  : weight = 1.07349
  elif(pt>=690  and pt < 740)  : weight = 1.10748
  elif(pt>=740  and pt < 790)  : weight = 1.06617
  elif(pt>=790  and pt < 840)  : weight = 1.05616
  elif(pt>=840  and pt < 900)  : weight = 1.1149
  elif(pt>=900  and pt < 960)  : weight = 1.03164
  elif(pt>=960  and pt < 1020) : weight = 1.06872
  elif(pt>=1020 and pt < 1090) : weight = 0.981645
  elif(pt>=1090 and pt < 1160) : weight = 0.81729
  elif(pt>=1160)              : weight = 0.924246

  return weight


def getRenUpZ(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 1.38864
  elif(pt>=170  and pt < 200)  : weight = 1.43444
  elif(pt>=200  and pt < 230)  : weight = 1.39055
  elif(pt>=230  and pt < 260)  : weight = 1.2934
  elif(pt>=260  and pt < 290)  : weight = 1.30618
  elif(pt>=290  and pt < 320)  : weight = 1.33453
  elif(pt>=320  and pt < 350)  : weight = 1.27239
  elif(pt>=350  and pt < 390)  : weight = 1.25981
  elif(pt>=390  and pt < 430)  : weight = 1.18738
  elif(pt>=430  and pt < 470)  : weight = 1.21686
  elif(pt>=470  and pt < 510)  : weight = 1.14351
  elif(pt>=510  and pt < 550)  : weight = 1.13305
  elif(pt>=550  and pt < 590)  : weight = 1.12346
  elif(pt>=590  and pt < 640)  : weight = 1.05153
  elif(pt>=640  and pt < 690)  : weight = 0.98056
  elif(pt>=690  and pt < 740)  : weight = 1.00651
  elif(pt>=740  and pt < 790)  : weight = 0.96642
  elif(pt>=790  and pt < 840)  : weight = 0.956161
  elif(pt>=840  and pt < 900)  : weight = 0.998079
  elif(pt>=900  and pt < 960)  : weight = 0.927421
  elif(pt>=960  and pt < 1020) : weight = 0.954094
  elif(pt>=1020 and pt < 1090) : weight = 0.883954
  elif(pt>=1090 and pt < 1160) : weight = 0.729348
  elif(pt>=1160)              : weight = 0.833531

  return weight


def getRenDownZ(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 1.55029
  elif(pt>=170  and pt < 200)  : weight = 1.64816
  elif(pt>=200  and pt < 230)  : weight = 1.59095
  elif(pt>=230  and pt < 260)  : weight = 1.48394
  elif(pt>=260  and pt < 290)  : weight = 1.49876
  elif(pt>=290  and pt < 320)  : weight = 1.55923
  elif(pt>=320  and pt < 350)  : weight = 1.49283
  elif(pt>=350  and pt < 390)  : weight = 1.49279
  elif(pt>=390  and pt < 430)  : weight = 1.39829
  elif(pt>=430  and pt < 470)  : weight = 1.46136
  elif(pt>=470  and pt < 510)  : weight = 1.38252
  elif(pt>=510  and pt < 550)  : weight = 1.36081
  elif(pt>=550  and pt < 590)  : weight = 1.37848
  elif(pt>=590  and pt < 640)  : weight = 1.28575
  elif(pt>=640  and pt < 690)  : weight = 1.17335
  elif(pt>=690  and pt < 740)  : weight = 1.21841
  elif(pt>=740  and pt < 790)  : weight = 1.17707
  elif(pt>=790  and pt < 840)  : weight = 1.16846
  elif(pt>=840  and pt < 900)  : weight = 1.25225
  elif(pt>=900  and pt < 960)  : weight = 1.15246
  elif(pt>=960  and pt < 1020) : weight = 1.20461
  elif(pt>=1020 and pt < 1090) : weight = 1.09389
  elif(pt>=1090 and pt < 1160) : weight = 0.921666
  elif(pt>=1160)              : weight = 1.02897

  return weight


def getFacUpZ(  pt):

  weight = 1.0
  if     (pt < 170)              : weight = 1.48202
  elif(pt>=170  and pt < 200)  : weight = 1.54551
  elif(pt>=200  and pt < 230)  : weight = 1.49889
  elif(pt>=230  and pt < 260)  : weight = 1.37781
  elif(pt>=260  and pt < 290)  : weight = 1.39565
  elif(pt>=290  and pt < 320)  : weight = 1.42918
  elif(pt>=320  and pt < 350)  : weight = 1.35817
  elif(pt>=350  and pt < 390)  : weight = 1.34966
  elif(pt>=390  and pt < 430)  : weight = 1.26377
  elif(pt>=430  and pt < 470)  : weight = 1.30245
  elif(pt>=470  and pt < 510)  : weight = 1.22426
  elif(pt>=510  and pt < 550)  : weight = 1.20589
  elif(pt>=550  and pt < 590)  : weight = 1.20726
  elif(pt>=590  and pt < 640)  : weight = 1.12468
  elif(pt>=640  and pt < 690)  : weight = 1.04086
  elif(pt>=690  and pt < 740)  : weight = 1.06898
  elif(pt>=740  and pt < 790)  : weight = 1.02577
  elif(pt>=790  and pt < 840)  : weight = 1.01591
  elif(pt>=840  and pt < 900)  : weight = 1.06486
  elif(pt>=900  and pt < 960)  : weight = 0.986133
  elif(pt>=960  and pt < 1020) : weight = 1.01224
  elif(pt>=1020 and pt < 1090) : weight = 0.937427
  elif(pt>=1090 and pt < 1160) : weight = 0.77774
  elif(pt>=1160)              : weight = 0.87435

  return weight


def getFacDownZ(  pt):
  weight = 1.0
  if     (pt < 170)              : weight = 1.47257
  elif(pt>=170  and pt < 200)  : weight = 1.54471
  elif(pt>=200  and pt < 230)  : weight = 1.49357
  elif(pt>=230  and pt < 260)  : weight = 1.41154
  elif(pt>=260  and pt < 290)  : weight = 1.42113
  elif(pt>=290  and pt < 320)  : weight = 1.46933
  elif(pt>=320  and pt < 350)  : weight = 1.41145
  elif(pt>=350  and pt < 390)  : weight = 1.40379
  elif(pt>=390  and pt < 430)  : weight = 1.32597
  elif(pt>=430  and pt < 470)  : weight = 1.37278
  elif(pt>=470  and pt < 510)  : weight = 1.2973
  elif(pt>=510  and pt < 550)  : weight = 1.28653
  elif(pt>=550  and pt < 590)  : weight = 1.28624
  elif(pt>=590  and pt < 640)  : weight = 1.20542
  elif(pt>=640  and pt < 690)  : weight = 1.11142
  elif(pt>=690  and pt < 740)  : weight = 1.1518
  elif(pt>=740  and pt < 790)  : weight = 1.11242
  elif(pt>=790  and pt < 840)  : weight = 1.1022
  elif(pt>=840  and pt < 900)  : weight = 1.17247
  elif(pt>=900  and pt < 960)  : weight = 1.08313
  elif(pt>=960  and pt < 1020) : weight = 1.13457
  elif(pt>=1020 and pt < 1090) : weight = 1.03169
  elif(pt>=1090 and pt < 1160) : weight = 0.862166
  elif(pt>=1160)              : weight = 0.981254
  return weight
