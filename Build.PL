use strict;
use Module::Build;

my $build = Module::Build->new(
	dist_name => 'DateTime-Util-Astro',
	dist_version => '0.01',
	license => 'perl',
	build_requires => {
		'Test::More' => 0
	},
	requires => {
		'DateTime' => 0,
		'DateTime::Util::Calc' => 0,
		'Params::Validate' => 0,
		'Math::Trig' => 0,
		'POSIX' => 0
	},
	recommends => {
		'perl' => '5.6.1',
		'Math::BigInt' => 0,
		'Math::BigFloat' => 0,
		'Math::BigInt::GMP' => 0,
	},
	sign => 1,
	create_readme => 0,
    create_makefile_pl => 'passthrough'
);
$build->create_build_script;
