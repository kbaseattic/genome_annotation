

# This is a simple script to test the annotation service.  The commands 
# are based on the tutorial 
# http://www.kbase.us/developer-zone/tutorials/getting-started/annotating-a-genome-using-kbase-tools/

#from subprocess import call
#from subprocess import Popen, PIPE
from os import system
import sys

TEST_DIR = "g.2860"

class Functional_Test():

    def __init__(self):
        pass

    def build_project(self):
        system('mkdir g.2860')


    def obtain_contigs(self):
        cmd = "echo 'kb|g.2860' | genomes_to_contigs | contigs_to_sequences > g.2860.contigs"
        system(cmd)
#        p1 = Popen(["echo", "\'kb|g.2860\'", '-c'], stdout=PIPE)
#        p2 = Popen(['genomes_to_contigs'], stdin=p1.stdout, stdout=PIPE)
#        p3 = Popen(['contigs_to_sequences', '>', 'g.2860/g.2860.contigs'], stdin=p2.stdout, stdout=PIPE)
#        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
#        output = p3.communicate()
#        print output

    def fasta_to_genome(self):
        cmd = "fasta_to_genome 'Geobacter sulfurreducens KN400' Bacteria 11 < g.2860.contigs \
                > genome"
        system(cmd)


    def annotate_genome(self):
        cmd = "annotate_genome < genome > annotated.genome"
        system(cmd)

    def gemome_to_features(self):
        cmd = "genomeTO_to_feature_data < annotated.genome > features.txt"
        system(cmd)

    def genome_to_reconstruction(self):
        cmd = "genomeTO_to_reconstructionTO < annotated.genome > reconstruction"
        system(cmd)

    def reconstruction_to_roles(self):
        cmd = "reconstructionTO_to_roles < reconstruction > roles"
        system(cmd)

    def reconstruction_to_subsystems(self):
        cmd = "reconstructionTO_to_subsystems < reconstruction > subsystems"
        system(cmd)

    def roles_used_in_models(self):
        cmd = " all_roles_used_in_models > all.roles"
        system(cmd)

    def a_and_b_roles(self):
        cmd = "a_and_b roles all.roles > roles.for.models"
        system(cmd)

    def missed_roles(self):
        cmd = "echo 'kb|g.9032' | genomes_to_fids CDS | fids_to_roles 2> /dev/null | cut -f 3 > \
                roles.in.g.9032"
        system(cmd)

    def roles_not_in_new_annotation(self):
        cmd = "a_and_b roles.in.g.9032 all.roles > roles.for.models.g.9032"
        system(cmd)
        cmd = "a_not_b roles.for.models roles.for.models.g.9032 > roles.to.search.for"
        system(cmd)        

    def tear_down_project(self):
        files = ["all.roles", "roles", "roles.to.search.for", "annotated.genome", 
                 "g.2860.contigs", "roles.for.models", "subsystems",
                 "features.txt", "genome", "roles.for.models.g.9032",
                 "reconstruction", "roles.in.g.9032"]  
        for file in files:
            system('rm '+file)




if __name__ == "__main__":
    test = Functional_Test()
    
    print "\nObtaining contigs..."
    test.obtain_contigs()
    print 'Done.'

    print "\nFasta to genome command..."
    test.fasta_to_genome()
    print 'Done.'

    print "\nAnnotate genome..."
    test.annotate_genome()
    print 'Done.'

    print "\nGenome to features file..."
    test.gemome_to_features() 
    print 'Done.'

    print "\nGenome to reconstruction..."
    test.genome_to_reconstruction()
    print 'Done.'

    print "\nReconstruction to roles..."
    test.reconstruction_to_roles()
    print 'Done.'

    print "\nReconstruction to subsystems..."
    test.reconstruction_to_subsystems()
    print 'Done.'

    print "\nRoles used in models..."
    test.roles_used_in_models()
    print 'Done.'

    print "\nIntersection of roles..."
    test.a_and_b_roles()
    print 'Done.'

    print "\nMissed roles..."
    test.missed_roles()
    print 'Done.'

    print "\nRoles not in new annotation..."
    test.roles_not_in_new_annotation()
    print 'Done.'

    print '\nTearing down poject...'
    test.tear_down_project()
    print 'Done.'
