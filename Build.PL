use strict;
use Module::Build;

my %requires = (
    'DateTime' => 0,
    'DateTime::Util::Calc' => '0.02',
    'Params::Validate' => 0,
    'Math::Trig' => 0,
    'POSIX' => 0
);

my $do_memoize = 0;
my $memoize_installed = eval { require Memoize };

if (!$memoize_installed) {
    $do_memoize = Module::Build->y_n(
      "  Would you like to install Memoize?\n\n" .
      "  DateTime::Event::Lunar by default tries to cache dates of new moons\n" .
      "  using Memoize (if available), which in some cases improves performance.\n\n" .
      "  Install Memoize? ",
      "n"
    );
}

if ($do_memoize) {
    print "Adding Memoize to list of prerequisites...\n";
    $requires{ Memoize } = 0;
}

my $build = Module::Build->new(
    dist_abstract => "Utilities for DateTime Extensions That Require Astronomical Calculations",
    dist_name => 'DateTime-Util-Astro',
    dist_version => '0.02',
    license => 'perl',
    build_requires => {
        'Test::More' => 0
    },
    requires => \%requires,
    recommends => {
        'perl' => '5.6.1',
    },
    sign => 1,
    create_makefile_pl => 'passthrough'
);
$build->create_build_script;
