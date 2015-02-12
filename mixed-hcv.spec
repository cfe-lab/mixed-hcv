# -*- mode: python -*-
from glob import glob

a = Analysis(['mixed-hcv.py'],
             pathex=['/Users/art/git/mixed-hcv'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)

#files = glob('data/*.bt2')
#for file in files:
#    print file
#    a.datas += [(file, file, 'DATA')]

pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='mixed-hcv',
          debug=False,
          strip=None,
          upx=True,
          console=True,
          icon='res/mixer.icns')
