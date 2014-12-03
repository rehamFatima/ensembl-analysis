# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
package Bio::EnsEMBL::Analysis::Tools::LayerFilter::AbstractLayerFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );



sub new{
  my ($class, @args) = @_;
  my $self = bless {},$class;

  if (scalar(@args)) {
    throw("$class should have no args in new");
  }

  return $self;
}

sub has_overlap {
    throw("You need to implement a method has_overlap!");
}

sub add_supporting_evidence {
    my ($gene, $discarded_gene) = @_;
}

sub filter {
  my ($self, $these, $others) = @_;

  # interference is judged by overlap at coding exon level
  # assumption is that @others is sorted by gene start

  my @filtered;

  my $cur_idx = 0;

  foreach my $obj (@$these) {
    my (@genomic_overlap, $left_bound);

    my $overlap = 0;
    for(my $i=$cur_idx; $i < @$others; $i++) {
        my $o_obj = $others->[$i];

        if ($o_obj->strand != $obj->strand or $o_obj->end < $obj->start) {
            next;
        } else {
            $left_bound = $i;
            if ($o_obj->start > $obj->end) {
                last;
            } else {
                $overlap = $self->has_overlap($obj, $o_obj);
            }
        }
    }

    $cur_idx = $left_bound if defined $left_bound;

    if (not $overlap) {
      push @filtered, $obj;
    }
  }

  return \@filtered;
}

1;
