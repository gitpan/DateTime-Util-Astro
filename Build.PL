# Build.PL,v 1.4 2005/01/07 12:18:59 lestrrat Exp
#
# Daisuke Maki <dmaki@cpan.org>
# All rights reserved.

use strict;
use Module::Build;

my $class = Module::Build->subclass(code => q|
    sub ACTION_manifest
    {
        my $self = shift;
        $self->SUPER::ACTION_manifest();
        if ($self->create_makefile_pl) {
            $self->_add_to_manifest('MANIFEST', 'Makefile.PL        Added here by Build.PL');
        }
    }
|);

my $build = $class->new(
    module_name => 'DateTime-Util-Astro',
    dist_version => '0.05',
    license => 'perl',
    build_requires => {
        'Test::More' => 0
    },
    requires => {
        'Cache::Cache'         => 0,
        'DateTime'             => 0,
        'DateTime::Util::Calc' => '0.08',
        'Math::BigFloat'       => 0,
        'Math::BigInt'         => 0,
        'Math::Trig'           => 0,
        'POSIX'                => 0
    },
    sign => 1,
    create_makefile_pl => 'traditional'
);
$build->create_build_script;
