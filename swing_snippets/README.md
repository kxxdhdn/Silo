# swing_snippets
For synthetic photometry

Installation Guide
-------------------

1. Open _utilities.f90_, search __/path/of/this/file__ and modify __LibF__ with the location of this file in your machine.
2. Launch a terminal window and change directories to the same directory.
3. Execute the installation program by entering the following command.
```bash
make
```
4. Add this directory to your PATH environment variable (e.g. in _.bashrc_)
5. [Optional] Run the test. (sss: Swing Snippets for Synthetic photometry)
```bash
cd Tests
python test_sss.py
```
6. To uninstall, just run 
```bash
make mrproper
```

Notices
--------

_swing_snippets_ is extracted from F. Galliano's _SwING_ library which is not yet published. 

Please report any bug or advices to dangning.hu@cea.fr

vlog
-----
- v0.1 (20200413)
  - Merged to astylo (see calib)
- v0 (20200310)
  - synphot_v0 (see /astylo/archives)