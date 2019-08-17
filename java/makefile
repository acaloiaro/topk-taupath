JAVAC = javac
N = 512
sources =  $(wildcard *.java) $(wildcard **/*.java)
classes = $(sources:.java=.class)
source_dir = ./tktp
manifest_file = $(source_dir)/MANIFEST.mf
jar_file = tktp.jar

all: compile jar clean_build
.PHONY : all

compile: $(classes)

%.class: %.java
	$(JAVAC) $<

jar:
	@echo "Manifest-Version: 1.0" > $(manifest_file)
	@echo "Class-Path: ." >> $(manifest_file)
	@echo "Main-Class: tktp.YuFastBCS" >> $(manifest_file)
	@echo "" >> $(manifest_file)

	jar -cmf $(manifest_file) $(jar_file) $(classes)

clean_build:
	rm -f $(source_dir)/*.class
	rm -f $(manifest_file)

clean_bin:
	rm -rf *.jar

clean: clean_build clean_bin

run:
	java -jar $(jar_file)

validate:
	$(eval X := $(shell Rscript --vanilla ../R/util_get_samples.r $(N)))
	$(eval Y := $(shell Rscript --vanilla ../R/util_get_samples.r $(N)))
	$(eval JAVAPI := $(shell java -cp tktp.jar tktp.ValidationHelper FastBCS "$(X)" "$(Y)"))
	$(eval RPI := $(shell Rscript --vanilla ../R/util_run_fastbcs.r "$(X)" "$(Y)" | xargs))

	@echo "$(JAVAPI)"
	@echo ""
	@echo "$(RPI)"

	@echo "validates?" $(shell [[ "$(RPI)" == "$(JAVAPI)" ]] && echo "yes" || echo "no")

validate_java:
	$(eval X := $(shell Rscript --vanilla ../R/util_get_samples.r $(N)))
	$(eval Y := $(shell Rscript --vanilla ../R/util_get_samples.r $(N)))
	$(eval YUPI := $(shell java -cp tktp.jar tktp.ValidationHelper FastBCS "$(X)" "$(Y)"))
	$(eval FBCS2PI := $(shell java -cp tktp.jar tktp.ValidationHelper FastBCS2 "$(X)" "$(Y)"))

	@echo $(YUPI)
	@echo ""
	@echo $(FBCS2PI)

	@echo "validates?" $(shell [[ "$(YUPI)" == "$(FBCS2PI)" ]] && echo "yes" || echo "no")
