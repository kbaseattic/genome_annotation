use Bio::KBase::CDMI::Client;
use strict;
use Data::Dumper;

my $c = Bio::KBase::CDMI::Client->new('http://kbase.us/services/cdmi_api');

my $g = 'kb|g.0';

my @roles;

my @xroles = (
'2-isopropylmalate synthase (EC 2.3.3.13)', 
	    '3-isopropylmalate dehydratase large subunit (EC 4.2.1.33)', 
	    '3-isopropylmalate dehydratase small subunit (EC 4.2.1.33)', 
	    '3-isopropylmalate dehydrogenase (EC 1.1.1.85)', 
	    );

if (!@roles)
{
    if (open(R, "<", "roles.tmp"))
    {
	@roles = <R>;
	chomp @roles;
	close(R);
    }
    else
    {
	my $fids = $c->genomes_to_fids([$g], ['CDS']);
	$fids = $fids->{$g};
	my $roles = $c->fids_to_roles($fids);
	my %roles;
	for my $v (values (%$roles))
	{
	    $roles{$_} = 1 foreach @$v;
	}
	@roles = sort keys %roles;
	open(R, ">", "roles.tmp") or die;
	print R "$_\n" foreach @roles;
	close(R);
    }
}

my $u = Bio::KBase::GenomeAnnotation::Utilities->new($c);

my $x = $u->metabolic_reconstruction(\@roles);
print join("\t", @$_), "\n" foreach @$x;

package Bio::KBase::GenomeAnnotation::Utilities;
use strict;
use Data::Dumper;

sub new
{
    my($class, $cdm) = @_;

    my $self = {
	cdm => $cdm,
    };

    return bless $self, $class;
    
}

sub metabolic_reconstruction {
    # Get the parameters.
    my ($self, $id_roles) = @_;

    my $retVal = [];

    # This counter will be used to generate user IDs for roles without them.
    my $next = 1000;

    my @id_roles1 = map { (ref $_ ? $_ : [$_, "FR" . ++$next]) } @$id_roles;

    my @id_roles = ();
    foreach my $tuple (@id_roles1)
    {
	my($function,$id) = @$tuple;
	foreach my $role (split(/(?:; )|(?: [\]\@] )/,$function))
	{
	    push(@id_roles,[$role,$id]);
	}
    }

    my %big;
    my $id_display = 1;
    map {push(@{$big{$_->[0]}}, $_->[1])} @id_roles;

    my $cdm = $self->{cdm};

    my %ss_roles;

    my $block = 2000;
    my $idx = 0;
    while (1)
    {
	my $slist = $cdm->all_entities_Subsystem($idx, $block, ['id']);
	last if %$slist == 0;
	my $ret = $cdm->get_relationship_Includes([keys %$slist], [], ['abbreviation', 'auxiliary', 'to_link'], []);
	for my $ent (@$ret)
	{
	    my($l_from, $l_rel, $l_to) = @$ent;
	    next if $l_rel->{auxiliary};
	    $ss_roles{$l_rel->{'from_link'}}->{$l_rel->{to_link}} = $l_rel->{abbreviation};
	}
	$idx += $block;
    }

    my @to_check;
    foreach my $sub (keys %ss_roles) {
	my $roles = $ss_roles{$sub};
	my @rolesubset = grep { $big{$_} } keys %$roles;
        my @abbr = map{$roles->{$_}} @rolesubset;
        my $set =  join(" ",  @abbr);
        if (@abbr > 0) {
	    push(@to_check, [$sub, $set]);
	}
    }
    
    #
    # Retrieve variants for selected subsystems.
    #

    # print Dumper([map { $_->[0]} @to_check]);
    
    my $res = $cdm->get_relationship_Describes([map { $_->[0]} @to_check], [], [], ['role_rule', 'code']);
    for my $ent (@$res)
    {
	my($f, $r, $t) = @$ent;
	# next unless $r->{from_link} eq '16S rRNA modification within P site of ribosome';
	my $rr = $t->{role_rule} || [];
	# print join("\t", $r->{from_link}, $t->{id}, $t->{code}, scalar @$rr), "\n";
    }

    my %ss_describes;
    for my $ent (@$res)
    {
	my($l_from, $l_rel, $l_to) = @$ent;
	push(@{$ss_describes{$l_rel->{from_link}}}, [@{$l_to}{'id', 'code', 'role_rule'}]);
    }

    #open(O, ">", "dump.$$");
    #print O Dumper($cdm, $res);
    #close(O);
    # print Dumper(\@to_check, \%ss_describes);

    for my $ent (@to_check)
    {
	my($sub, $set) = @$ent;


	my $roles = $ss_roles{$sub};
	my ($variant, $size) = $self->get_max_subset($sub, $set, $ss_describes{$sub});
	if ($variant) {
	    foreach my $role (keys %$roles) {
		if ($id_display) {
		    if (exists $big{$role}) {
			foreach my $id (@{$big{$role}}) {
			    push (@$retVal, [$sub, $variant, $role, $id]);
			}
		    }
		} else {
		    push (@$retVal, [$sub, $variant, $role]);
		}
	    }
	}
    }
    # Return the result.
    return $retVal;
}

=head3 get_max_subset

    my ($ss, $max_variant, $max_size) = $ssObject->get_max_subset($sub, $setA, $variants);

Given a subsystem ID and a role rule, return the ID of the variant for
the subsystem that matches the most roles in the rule and the number of
roles matched.

=over 4

=item sub

Name (ID) of the subsystem whose variants are to be examined.

=item setA

A space-delimited list of role abbreviations, lexically ordered. This provides
a unique specification of the roles in the set.

=item variants

List of triples [id, code, code-rule] describing the variants for this subsystem.

=item RETURN

Returns a 2-element list consisting of name variant found (subsystem name, colon,
and variant code) and the number of roles matched.

=back

=cut

sub get_max_subset {
    my ($self, $sub, $setA, $variants) = @_;
    my $max_size = 0;
    my $max_set;
    my $max_variant;
    my %set_hash;

    for my $v (@$variants)
    {
	my($id, $variantCode, $variantRoleRuleList) = @$v;

        foreach my $setB (@$variantRoleRuleList) {
	    my $size = is_A_a_superset_of_B($setA, $setB);
	    # printf "%3d: setA=$setA\n     setB=$setB\n", $size;
	    if ($size  && $size > $max_size) {
		$max_size = $size;
		$max_set = $setB;
		$max_variant = $variantCode;
	    }
        }
    }
#    if ($max_size) {
#	print STDERR "Success $max_variant $max_set $max_size\n";
#    }
    return ($max_variant, $max_size);
}


=head3 is_A_a_superset_of_B

    my $size = SS::is_A_a_superset_of_B($a, $b);

This method takes as input two role rules, and returns 0 if the first
role rule is NOT a superset of the second; otherwise, it returns the size
of the second rule. A role rule is a space-delimited list of role
abbreviations in lexical order. This provides a unique identifier for a
set of roles in a subsystem.

=over 4

=item a

First role rule.

=item b

Second role rule.

=item RETURN

Returns 0 if the first rule is NOT a superset of the second and the size of the
second rule if it is. As a result, if the first rule IS a superset, this method
will evaluate to TRUE, and to FALSE otherwise.

=back

=cut

sub is_A_a_superset_of_B {
    my ($a, $b) = @_;
    my @a = split(" ", $a);
    my @b = split(" ", $b);
    if (@b > @a) {
            return(0);
    }
    my %given;
    map { $given{$_} = 1} @a;
    map { if (! $given{$_}) {return 0}} split(" ", $b);
    my $l = scalar(@b);
    return scalar(@b);
}


1;
