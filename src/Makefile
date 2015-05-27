install: sbams_mvlr/sbams_mvlr  sbams_sslr/sbams_sslr
	mkdir -p fmqtl_bin
	cp sbams_mvlr/sbams_mvlr fmqtl_bin/
	cp sbams_sslr/sbams_sslr fmqtl_bin/
	cp utility/*.pl utility/fm_qtl fmqtl_bin/
	make clean
sbams_mvlr/sbams_mvlr:
	cd sbams_mvlr;make;cd ..
sbams_sslr/sbams_sslr:
	cd sbams_sslr;make ;cd ..
clean:
	cd sbams_mvlr;make clean;cd ..
	cd sbams_sslr;make clean;cd ..