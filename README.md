### Dependencies

* External Package

* Ports that need to be open
  Port 7050

### Setup using the kbase VMs
=======
0.  Start the VM and clone the git repo.
    nova boot .... (options will change over time)
    ssh ubuntu@<vm host>

1. Following an updated version of the directions
from: https://trac.kbase.us/projects/kbase/wiki/IntegrationTargets
   sudo bash
   cd /kb
   git clone kbase@git.kbase.us:/dev_container.git
   cd dev_container/modules
   git clone kbase@git.kbase.us:/genome_annotation.git
   cd ..
   ./bootstrap /kb/runtime
   . user-env.sh

2. The make target deploy will install the service, libraries and documentation
   cd modules/genome_annotation
   make deploy

3. Run the internal unit tests
   make test   

4. The service is started/stopped by using the start_service and stop_service
scripts in /kb/deployment/services/genome_annotation.
   cd /kb/deployment/services/genome_annotation
   ./start_service # or ./stop_service
