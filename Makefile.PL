#!perl
use strict;
use inc::Module::Install;

name('DateTime-Util-Astro');
license('perl');
version_from('lib/DateTime/Util/Astro.pm');
requires 'DateTime';
requires 'DateTime::Util::Calc' => '0.12';
requires 'Math::BigInt';
requires 'Math::Trig';
requires 'POSIX';

auto_install;
WriteAll;