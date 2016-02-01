TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime
SERVER_SPEC = GenomeAnnotation.spec

SERVICE_MODULE = lib/Bio/KBase/GenomeAnnotation/Service.pm

SERVICE = genome_annotation
SERVICE_PORT = 7050
SERVICE_ALT_PORT = 7136

SERVICE_URL = https://kbase.us/services/$(SERVICE)

SERVICE_NAME = GenomeAnnotation
SERVICE_NAME_PY = $(SERVICE_NAME)

SERVICE_PSGI_FILE = $(SERVICE_NAME).psgi

SRC_SERVICE_PERL = $(wildcard service-scripts/*.pl)
BIN_SERVICE_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PERL))))
DEPLOY_SERVICE_PERL = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PERL))))


ifdef TEMPDIR
TPAGE_TEMPDIR = --define kb_tempdir=$(TEMPDIR)
endif

ifdef SERVICE_ALT_PORT
TPAGE_SERVICE_ALT_PORT = --define kb_service_alt_port=$(SERVICE_ALT_PORT) 
endif

TPAGE_ARGS = --define kb_top=$(TARGET) \
	--define kb_runtime=$(DEPLOY_RUNTIME) \
	--define kb_service_name=$(SERVICE) \
	--define kb_service_port=$(SERVICE_PORT) \
	--define kb_service_port=$(SERVICE_ALT_PORT) \
	$(TPAGE_SERVICE_ALT_PORT) \
	$(TPAGE_TEMPDIR)

TESTS = $(wildcard t/client-tests/*.t)

all: bin compile-typespec service

test:
	# run each test
	echo "RUNTIME=$(DEPLOY_RUNTIME)\n"
	for t in $(TESTS) ; do \
		if [ -f $$t ] ; then \
			$(DEPLOY_RUNTIME)/bin/perl $$t ; \
			if [ $$? -ne 0 ] ; then \
				exit 1 ; \
			fi \
		fi \
	done

service: $(SERVICE_MODULE)

$(SERVICE_MODULE): $(SERVER_SPEC) 
	./recompile_typespec 

compile-typespec: Makefile
	mkdir -p lib/biokbase/$(SERVICE_NAME_PY)
	touch lib/biokbase/__init__.py #do not include code in biokbase/__init__.py
	touch lib/biokbase/$(SERVICE_NAME_PY)/__init__.py 
	mkdir -p lib/javascript/$(SERVICE_NAME)
	compile_typespec \
		--psgi $(SERVICE_PSGI_FILE) \
		--impl Bio::KBase::$(SERVICE_NAME)::%sImpl \
		--service Bio::KBase::$(SERVICE_NAME)::Service \
		--client Bio::KBase::$(SERVICE_NAME)::Client \
		--py biokbase/$(SERVICE_NAME_PY)/client \
		--js javascript/$(SERVICE_NAME)/Client \
		--url $(SERVICE_URL) \
		--enable-retries \
		$(SERVER_SPEC) lib
	-rm -f lib/$(SERVER_MODULE)Server.py
	-rm -f lib/$(SERVER_MODULE)Impl.py
	-rm -f lib/CDMI_EntityAPIImpl.py

bin: $(BIN_PERL) $(BIN_DIR)/kmer_guts $(BIN_SERVICE_PERL)

$(BIN_DIR)/kmer_guts: src/kmer_guts
	rm -f $(BIN_DIR)/kmer_guts
	cp src/kmer_guts $(BIN_DIR)/kmer_guts

src/kmer_guts: src/kmer_guts.c
	cd src; $(CC) $(CFLAGS) -O -o kmer_guts kmer_guts.c

deploy: deploy-client deploy-service
deploy-all: deploy-client deploy-service
deploy-client: compile-typespec deploy-docs deploy-libs deploy-scripts deploy-alt-scripts


RAST_SCRIPTS = $(filter scripts/rast2-%,$(SRC_PERL))
deploy-alt-scripts:
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PERL_PATH=$(TARGET)/lib ; \
	for src in $(RAST_SCRIPTS) ; do \
		basefile=`basename $$src`; \
		base=`basename $$src .pl`; \
		alt=`echo $$base | sed 's/^rast2-/rast-/'`; \
		echo install $$src $$alt ; \
		cp $$src $(TARGET)/plbin ; \
		$(WRAP_PERL_SCRIPT) "$(TARGET)/plbin/$$basefile" $(TARGET)/bin/$$alt ; \
	done 

deploy-guts: deploy-dir
	rm -f $(TARGET)/services/$(SERVICE)/bin/kmer_guts
	cp $(BIN_DIR)/kmer_guts $(TARGET)/services/$(SERVICE)/bin/kmer_guts

deploy-service: deploy-dir deploy-monit deploy-libs deploy-guts deploy-service-scripts
	$(TPAGE) $(TPAGE_ARGS) service/start_service.tt > $(TARGET)/services/$(SERVICE)/start_service
	chmod +x $(TARGET)/services/$(SERVICE)/start_service
	$(TPAGE) $(TPAGE_ARGS) service/stop_service.tt > $(TARGET)/services/$(SERVICE)/stop_service
	chmod +x $(TARGET)/services/$(SERVICE)/stop_service

deploy-service-scripts:
	export KB_TOP=$(TARGET); \
	export KB_RUNTIME=$(DEPLOY_RUNTIME); \
	export KB_PERL_PATH=$(TARGET)/lib ; \
	export PATH_PREFIX=$(TARGET)/services/$(SERVICE)/bin:$(TARGET)/services/cdmi_api/bin; \
	for src in $(SRC_SERVICE_PERL) ; do \
	        basefile=`basename $$src`; \
	        base=`basename $$src .pl`; \
	        echo install $$src $$base ; \
	        cp $$src $(TARGET)/plbin ; \
	        $(WRAP_PERL_SCRIPT) "$(TARGET)/plbin/$$basefile" $(TARGET)/services/$(SERVICE)/bin/$$base ; \
	done

deploy-monit:
	$(TPAGE) $(TPAGE_ARGS) service/process.$(SERVICE).tt > $(TARGET)/services/$(SERVICE)/process.$(SERVICE)

deploy-docs:
	-mkdir doc
	-mkdir $(SERVICE_DIR)
	-mkdir $(SERVICE_DIR)/webroot
	mkdir -p doc
	$(DEPLOY_RUNTIME)/bin/pod2html -t "Genome Annotation Service API" lib/Bio/KBase/GenomeAnnotation/GenomeAnnotationImpl.pm > doc/genomeanno_impl.html
	cp doc/*html $(SERVICE_DIR)/webroot/.

deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/webroot ] ; then mkdir $(SERVICE_DIR)/webroot ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi

$(BIN_DIR)/%: service-scripts/%.pl $(TOP_DIR)/user-env.sh
	$(WRAP_PERL_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

$(BIN_DIR)/%: service-scripts/%.py $(TOP_DIR)/user-env.sh
	$(WRAP_PYTHON_SCRIPT) '$$KB_TOP/modules/$(CURRENT_DIR)/$<' $@

include $(TOP_DIR)/tools/Makefile.common.rules
