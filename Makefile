TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

TARGET ?= /kb/deployment
DEPLOY_RUNTIME ?= /kb/runtime
SERVER_SPEC = GenomeAnnotation.spec

SERVICE_MODULE = lib/Bio/KBase/GenomeAnnotation/Service.pm

SERVICE = genome_annotation
SERVICE_PORT = 7050

TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(DEPLOY_RUNTIME) --define kb_service_name=$(SERVICE) \
	--define kb_service_port=$(SERVICE_PORT)

TESTS = $(wildcard t/client-tests/*.t)

all: bin service

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

bin: $(BIN_PERL)

deploy: deploy-client deploy-service
deploy-all: deploy-client deploy-service
deploy-client: deploy-docs deploy-libs deploy-scripts

deploy-libs: recompile-typespec

recompile-typespec:
	./recompile_typespec

deploy-service: deploy-monit deploy-libs
	$(TPAGE) $(TPAGE_ARGS) service/start_service.tt > $(TARGET)/services/$(SERVICE)/start_service
	chmod +x $(TARGET)/services/$(SERVICE)/start_service
	$(TPAGE) $(TPAGE_ARGS) service/stop_service.tt > $(TARGET)/services/$(SERVICE)/stop_service
	chmod +x $(TARGET)/services/$(SERVICE)/stop_service

deploy-monit:
	$(TPAGE) $(TPAGE_ARGS) service/process.$(SERVICE).tt > $(TARGET)/services/$(SERVICE)/process.$(SERVICE)

deploy-docs:
	-mkdir doc
	-mkdir $(SERVICE_DIR)
	-mkdir $(SERVICE_DIR)/webroot
	mkdir -p doc
	$(DEPLOY_RUNTIME)/bin/pod2html -t "Genome Annotation Service API" lib/Bio/KBase/GenomeAnnotation/GenomeAnnotationImpl.pm > doc/genomeanno_impl.html
	cp doc/*html $(SERVICE_DIR)/webroot/.

include $(TOP_DIR)/tools/Makefile.common.rules
