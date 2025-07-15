TOPDIR=..

include $(TOPDIR)/build/Makefile.inc


default: test coef_io other mw_scatt mw_scatt_coef coef_io_11 emis_atlas brdf_atlas wrapper gui 


all: bin

bin: main/bin parallel/bin hdf/bin coef_io/bin mw_scatt/bin other/bin coef_io_11/bin brdf_atlas/bin emis_atlas/bin mw_scatt_coef/bin gui/bin test/bin wrapper/bin


lib: main/lib parallel/lib hdf/lib coef_io/lib mw_scatt/lib other/lib coef_io_11/lib brdf_atlas/lib emis_atlas/lib mw_scatt_coef/lib gui/lib test/lib wrapper/lib


intf: main/intf parallel/intf hdf/intf coef_io/intf mw_scatt/intf other/intf coef_io_11/intf brdf_atlas/intf emis_atlas/intf mw_scatt_coef/intf gui/intf test/intf wrapper/intf


dir: main/dir parallel/dir hdf/dir coef_io/dir mw_scatt/dir other/dir coef_io_11/dir brdf_atlas/dir emis_atlas/dir mw_scatt_coef/dir gui/dir test/dir wrapper/dir



mw_scatt: mw_scatt/bin

mw_scatt/bin: mw_scatt/lib
	cd mw_scatt && $(MAKE) bin

mw_scatt/lib: mw_scatt/intf hdf/lib coef_io/lib main/lib
	cd mw_scatt && $(MAKE) lib

mw_scatt/intf: mw_scatt/dir
	cd mw_scatt && $(MAKE) intf

mw_scatt/dir:
	cd mw_scatt && $(MAKE) dir

mw_scatt/clean:
	cd mw_scatt && $(MAKE) clean

main: main/bin

main/bin: main/lib
	cd main && $(MAKE) bin

main/lib: main/intf 
	cd main && $(MAKE) lib

main/intf: main/dir
	cd main && $(MAKE) intf

main/dir:
	cd main && $(MAKE) dir

main/clean:
	cd main && $(MAKE) clean

gui: gui/bin

gui/bin: gui/lib
	cd gui && $(MAKE) bin

gui/lib: gui/intf emis_atlas/lib hdf/lib parallel/lib other/lib brdf_atlas/lib coef_io/lib main/lib
	cd gui && $(MAKE) lib

gui/intf: gui/dir
	cd gui && $(MAKE) intf

gui/dir:
	cd gui && $(MAKE) dir

gui/clean:
	cd gui && $(MAKE) clean

coef_io: coef_io/bin

coef_io/bin: coef_io/lib
	cd coef_io && $(MAKE) bin

coef_io/lib: coef_io/intf main/lib hdf/lib
	cd coef_io && $(MAKE) lib

coef_io/intf: coef_io/dir
	cd coef_io && $(MAKE) intf

coef_io/dir:
	cd coef_io && $(MAKE) dir

coef_io/clean:
	cd coef_io && $(MAKE) clean

other: other/bin

other/bin: other/lib
	cd other && $(MAKE) bin

other/lib: other/intf coef_io/lib main/lib hdf/lib
	cd other && $(MAKE) lib

other/intf: other/dir
	cd other && $(MAKE) intf

other/dir:
	cd other && $(MAKE) dir

other/clean:
	cd other && $(MAKE) clean

coef_io_11: coef_io_11/bin

coef_io_11/bin: coef_io_11/lib
	cd coef_io_11 && $(MAKE) bin

coef_io_11/lib: coef_io_11/intf main/lib coef_io/lib hdf/lib
	cd coef_io_11 && $(MAKE) lib

coef_io_11/intf: coef_io_11/dir
	cd coef_io_11 && $(MAKE) intf

coef_io_11/dir:
	cd coef_io_11 && $(MAKE) dir

coef_io_11/clean:
	cd coef_io_11 && $(MAKE) clean

mw_scatt_coef: mw_scatt_coef/bin

mw_scatt_coef/bin: mw_scatt_coef/lib
	cd mw_scatt_coef && $(MAKE) bin

mw_scatt_coef/lib: mw_scatt_coef/intf hdf/lib mw_scatt/lib coef_io/lib main/lib
	cd mw_scatt_coef && $(MAKE) lib

mw_scatt_coef/intf: mw_scatt_coef/dir
	cd mw_scatt_coef && $(MAKE) intf

mw_scatt_coef/dir:
	cd mw_scatt_coef && $(MAKE) dir

mw_scatt_coef/clean:
	cd mw_scatt_coef && $(MAKE) clean

test: test/bin

test/bin: test/lib
	cd test && $(MAKE) bin

test/lib: test/intf parallel/lib hdf/lib emis_atlas/lib mw_scatt/lib coef_io/lib main/lib other/lib brdf_atlas/lib
	cd test && $(MAKE) lib

test/intf: test/dir
	cd test && $(MAKE) intf

test/dir:
	cd test && $(MAKE) dir

test/clean:
	cd test && $(MAKE) clean

brdf_atlas: brdf_atlas/bin

brdf_atlas/bin: brdf_atlas/lib
	cd brdf_atlas && $(MAKE) bin

brdf_atlas/lib: brdf_atlas/intf coef_io/lib main/lib hdf/lib
	cd brdf_atlas && $(MAKE) lib

brdf_atlas/intf: brdf_atlas/dir
	cd brdf_atlas && $(MAKE) intf

brdf_atlas/dir:
	cd brdf_atlas && $(MAKE) dir

brdf_atlas/clean:
	cd brdf_atlas && $(MAKE) clean

parallel: parallel/bin

parallel/bin: parallel/lib
	cd parallel && $(MAKE) bin

parallel/lib: parallel/intf main/lib
	cd parallel && $(MAKE) lib

parallel/intf: parallel/dir
	cd parallel && $(MAKE) intf

parallel/dir:
	cd parallel && $(MAKE) dir

parallel/clean:
	cd parallel && $(MAKE) clean

wrapper: wrapper/bin

wrapper/bin: wrapper/lib
	cd wrapper && $(MAKE) bin

wrapper/lib: wrapper/intf mw_scatt/lib coef_io/lib main/lib brdf_atlas/lib other/lib parallel/lib hdf/lib emis_atlas/lib
	cd wrapper && $(MAKE) lib

wrapper/intf: wrapper/dir
	cd wrapper && $(MAKE) intf

wrapper/dir:
	cd wrapper && $(MAKE) dir

wrapper/clean:
	cd wrapper && $(MAKE) clean

hdf: hdf/bin

hdf/bin: hdf/lib
	cd hdf && $(MAKE) bin

hdf/lib: hdf/intf main/lib
	cd hdf && $(MAKE) lib

hdf/intf: hdf/dir
	cd hdf && $(MAKE) intf

hdf/dir:
	cd hdf && $(MAKE) dir

hdf/clean:
	cd hdf && $(MAKE) clean

emis_atlas: emis_atlas/bin

emis_atlas/bin: emis_atlas/lib
	cd emis_atlas && $(MAKE) bin

emis_atlas/lib: emis_atlas/intf main/lib coef_io/lib hdf/lib
	cd emis_atlas && $(MAKE) lib

emis_atlas/intf: emis_atlas/dir
	cd emis_atlas && $(MAKE) intf

emis_atlas/dir:
	cd emis_atlas && $(MAKE) dir

emis_atlas/clean:
	cd emis_atlas && $(MAKE) clean

clean:
	\rm -f $(OBJDIR)/*$(EXT_OBJ) $(LIBDIR)/*$(EXT_LIB) $(LIBDIR)/*$(EXT_SO) $(MODDIR)/*$(EXT_MOD) $(BINDIR)/*$(EXT_BIN) $(INCDIR)/*
	\rm -rf $(TMPDIR)/*

dist:
	cd .. ; tar cf - ReleaseNote.pdf docs gui src build data brdf_data emis_data rtcoef_rttov13 rttov_test wrapper --exclude src/dev --exclude src/lbl --exclude src/lbl-add --exclude src/baran --exclude src/coef_distribute --exclude "src/*/*.pl" --exclude "src/*/parallel.sh" --exclude "src/*/rttov_parallel*.txt" --exclude "docs/*lbl*" --exclude "docs/rttov-compile.*" --exclude "data/opac_interp_rh.py" --exclude "rtcoef_rttov13/extract*srf*.*" --exclude "rtcoef_rttov13/plot_lbl_*.py" --exclude "docs/*.odt" --exclude "docs/*.doc*" --exclude "docs/*.bib" --exclude "docs/*.tex" --exclude "docs/rttov_developer_guide.*" --exclude "docs/doxygen_doc_*" --exclude "docs/sharp_1983_*" --exclude gui/rttov_tests.tar.bz2 --exclude rttov_test/rttov_remote.pl --exclude "rttov_test/*_dev*" --exclude "rttov_test/*.1.*" --exclude "rttov_test/*.tar.gz" --exclude "rtcoef_rttov13/*/*.H5" --exclude "rtcoef_rttov13/rttov*pred*/*airs*" --exclude "rtcoef_rttov13/rttov*pred*/*iasi*" --exclude "rtcoef_rttov13/rttov*pred*/*cris*" --exclude "rtcoef_rttov13/rttov*pred*/*iris*" --exclude "rtcoef_rttov13/rttov*pred*/*mtg*irs*" --exclude "rtcoef_rttov13/cldaer*/*" --exclude "rtcoef_rttov13/mfasis_lut/*" --exclude "rtcoef_rttov13/hydrotable*/*.dat" --exclude "rtcoef_rttov13/pc/*.dat"  --exclude "rtcoef_rttov13/htfrtc/*.dat" --exclude "emis_data/*" --exclude "brdf_data/*" --exclude "wrapper/doxygen_doc_*" --exclude "wrapper/*_dev*" | xz > rttov_dist.tar.xz

distdev:
	cd .. ; tar cf - ReleaseNote.pdf docs gui src build data brdf_data emis_data rttov_test wrapper --exclude "docs/doxygen_doc_*" --exclude "rttov_test/*.tar.gz" --exclude "emis_data/*" --exclude "brdf_data/*" --exclude "wrapper/doxygen_doc_*" | gzip > rttov_dist_dev.tgz

distlbl:
	cd .. ; tar cf - docs/*lbl* src/lbl src/lbl-add rttov_lbl-data --exclude rttov_lbl-data/scripts_pascal --exclude rttov_lbl-data/rttov_lbl_instrument.kcarta --exclude rttov_lbl-data/instrument-data.kcarta | gzip > rttov_dist_lbl.tgz
