# This is a simple script to test the annotation service.  The commands 
# are based on the tutorial 
# http://www.kbase.us/developer-zone/tutorials/getting-started/annotating-a-genome-using-kbase-tools/

from os import system
import json
import sys

def check_columns(file_name, num_of_col):
    for line in open(file_name):
        row_entities = line.split('\t')
        assert len(row_entities) == num_of_col

def verify_json_file(file_name):
    json_file = open(file_name)
    try:
        data = json.load(json_file)
    except:
        print >> sys.stderr, file_name, "is not in json format!"


class FunctionalTest():
    def __init__(self):
        pass

    def build_project(self):
        system('mkdir g.2860')

    def obtain_contigs(self):
        cmd = "echo 'kb|g.2860' | genomes_to_contigs | contigs_to_sequences > g.2860.contigs"
        system(cmd)

    def fasta_to_genome(self):
        cmd = "fasta_to_genome 'Geobacter sulfurreducens KN400' Bacteria 11 < g.2860.contigs \
                > genome"
        system(cmd)

    def annotate_genome(self):
        cmd = "annotate_genome < genome > annotated.genome"
        system(cmd)

    def genome_to_features(self):
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
    test = FunctionalTest()
    
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
    test.genome_to_features() 
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
    print '***'

    # Verify files
    print "\nChecking columns features.txt..."
    check_columns('features.txt', 6)
    print 'Done.'

    print "\nVerifying reconstruction json..."
    verify_json_file('reconstruction')
    print 'Done.'

    print "\nChecking columns for roles..."
    check_columns('roles', 1)
    print 'Done.'

    print "\nChecking columns for subsytems..."
    check_columns('subsystems', 2)
    print 'Done.'

    print "\nChecking columns for roles.in.g.9032..."
    check_columns('roles.in.g.9032', 1)
    print 'Done.'

    print "\nChecking columns for subsytems..."
    check_columns('roles.to.search.for', 1)
    print 'Done.'


    print '\nTearing down poject...'
    test.tear_down_project()
    print 'Done.'
