#!perl
use strict;
use inc::Module::Install;

name('DateTime-Util-Astro');
license('perl');
version('0.08');
requires 'DateTime';
requires 'DateTime::Util::Calc' => '0.12';
requires 'Math::BigInt';
requires 'Math::Trig';
requires 'POSIX';

auto_install;
WriteAll;