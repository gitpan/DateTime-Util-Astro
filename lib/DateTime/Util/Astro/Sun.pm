package DateTime::Util::Astro::Sun;
use strict;
use vars qw($VERSION @ISA @EXPORT_OK);
BEGIN
{
	$VERSION = '0.01';
    @ISA = qw(Exporter);
    @EXPORT_OK = qw(
        solar_longitude
        solar_longitude_after
        solar_longitude_before
        estimate_prior_solar_longitude
    );
}

use DateTime;
use DateTime::Util::Calc qw(
    polynomial mod min max bf_downgrade angle bigfloat binary_search
    moment dt_from_moment sin_deg cos_deg tan_deg asin_deg acos_deg);
use DateTime::Util::Astro::Common
	qw(julian_centuries aberration nutation MEAN_TROPICAL_YEAR);
use POSIX();


# [1] Table 12.1 p.183 (zero-padded to align veritcally)
use constant SOLAR_LONGITUDE_ARGS => [
    #      X          Y             Z
    # left side of table 12.1
    [ 403406, 270.54861,      0.9287892 ],
    [ 119433,  63.91854,  35999.4089666 ],
    [   3891, 317.84300,  71998.2026100 ],
    [   1721, 240.05200,  36000.3572600 ],
    [    350, 247.23000,  32964.4678000 ],
    [    314, 297.82000, 445267.1117000 ],
    [    242, 166.79000,      3.1008000 ],
    [    158,   3.50000,    -19.9739000 ],
    [    129, 182.95000,   9038.0293000 ],
    [     99,  29.80000,  33718.1480000 ],
    [     86, 249.20000,  -2280.7730000 ],
    [     72, 257.80000,  31556.4930000 ],
    [     64,  69.90000,   9037.7500000 ],
    [     38, 197.10000,  -4444.1760000 ],
    [     32,  65.30000,  67555.3160000 ],
    [     28, 341.50000,  -4561.5400000 ],
    [     27,  98.50000,   1221.6550000 ],
    [     24, 110.00000,  31437.3690000 ],
    [     21, 342.60000, -31931.7570000 ],
    [     18, 256.10000,   1221.9990000 ],
    [     14, 242.90000,  -4442.0390000 ],
    [     13, 151.80000,    119.0660000 ],
    [     12,  53.30000,     -4.5780000 ],
    [     10, 205.70000,    -39.1270000 ],
    [     10, 146.10000,  90073.7780000 ],

    # right side of table 12.1
    [ 195207, 340.19128,  35999.1376958 ],
    [ 112392, 331.26220,  35998.7287385 ],
    [   2819,  86.63100,  71998.4403000 ],
    [    660, 310.26000,  71997.4812000 ],
    [    334, 260.87000,    -19.4410000 ],
    [    268, 343.14000,  45036.8840000 ],
    [    234,  81.53000,  22518.4434000 ],
    [    132, 132.75000,  65928.9345000 ],
    [    114, 162.03000,   3034.7684000 ],
    [     93, 266.40000,   3034.4480000 ],
    [     78, 157.60000,  29929.9920000 ],
    [     68, 185.10000,    149.5880000 ],
    [     46,   8.00000, 107997.4050000 ],
    [     37, 250.40000,    151.7710000 ],
    [     29, 162.70000,  31556.0800000 ],
    [     27, 291.60000, 107996.7060000 ],
    [     25, 146.70000,  62894.1670000 ],
    [     21,   5.20000,  14578.2980000 ],
    [     20, 230.90000,  34777.2430000 ],
    [     17,  45.30000,  62894.5110000 ],
    [     13, 115.20000, 107997.9090000 ],
    [     13, 285.30000,  16859.0710000 ],
    [     10, 126.60000,  26895.2920000 ],
    [     10,  85.90000,  12297.5360000 ]
];

# [1] p.182
sub solar_longitude
{
    my($dt) = Params::Validate::validate_pos(@_, { isa => 'DateTime' });

    my $c = julian_centuries($dt);
    my $big_ugly_number = 0;
    foreach my $numbers (@{ SOLAR_LONGITUDE_ARGS() }) {
        $big_ugly_number += bigfloat($numbers->[0]) * 
            sin_deg($numbers->[1] + $numbers->[2] * $c)
    }

    my $longitude = 282.7771834 +
        36000.76953744 * $c + 
        0.000005729577951308232 * $big_ugly_number;

    return mod($longitude + aberration($dt) + nutation($dt), 360);
}

# [1] p.184
sub solar_longitude_after
{
    my($dt, $phi) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { type => Params::Validate::SCALAR() }, 
    );

    my $epsilon = 10 ** -5;
    my $rate    = MEAN_TROPICAL_YEAR / 360;
    my $tau     = moment($dt) +
        $rate * mod($phi - solar_longitude($dt), 360);
    my $l       = max(moment($dt), $tau - 5);
    my $u       = $tau + 5;

print "approx date: ", dt_from_moment($tau)->datetime, "\n";
    my $rv = binary_search($l, $u,
        sub { abs($_[0] - $_[1]) <= $epsilon },
        sub { mod(solar_longitude(
            dt_from_moment($_[0])) - $phi, 360) < 180 } );
    return dt_from_moment(bf_downgrade($rv));
}

sub solar_longitude_before
{
    my($dt, $phi) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { type => Params::Validate::SCALAR() }, 
    );

    my $epsilon = 10 ** -5;
    my $rate    = MEAN_TROPICAL_YEAR / 360;
    my $tau     = moment($dt) +
        $rate * mod(solar_longitude($dt) - $phi, 360);
    my $l       = $tau - 5;
    my $u       = min(moment($dt), $tau + 5);


    my $rv = binary_search($l, $u,
        sub { abs($_[0] - $_[1]) <= $epsilon },
        sub { mod(solar_longitude(
            dt_from_moment($_[0])), 360) < 180 } );
    return dt_from_moment(bf_downgrade($rv));
}

# [1] p.203
sub estimate_prior_solar_longitude
{
    my($dt, $phi) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { type => Params::Validate::SCALAR() }, 
    );

    my $rate = MEAN_TROPICAL_YEAR / 360;
    my $tau  = moment($dt) - $rate *
        mod(solar_longitude($dt) - $phi, 360);
    my $delta = mod(solar_longitude(
        dt_from_moment($tau)) - $phi + 180, 360) - 180;

    my $rv = min(moment($dt), $tau - $rate * $delta);

    return dt_from_moment(bf_downgrade($rv));
    
}

# [1] p.198, errata 179
sub sine_offset
{
    my($dt, $location, $alpha) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { isa => 'Astro::Earth::Location' }, 
        { type => Params::Validate::SCALAR() } );
    my $phi   = $location->latitude;
    my $delta = asin_deg(
        sin_deg(obliquity($dt)) *
        sin_deg(solar_longitude($dt))
    );
    return tan_deg($phi) * tan_deg($delta) +
        (sin_deg($alpha) / (cos_deg($delta) * cos_deg($phi)));
}

# [1] p.198, errata 179
sub approx_moment_of_depression
{
    my($dt, $location, $alpha, $morning) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { isa => 'Astro::Earth::Location' }, 
        { type => Params::Validate::SCALAR() },
        { type => Params::Validate::BOOLEAN() } );

    my $try  = sine_offset($dt, $location, $alpha);
    my $date = POSIX::floor(moment($dt));
    my $alt  = $alpha >= 0 ?
        ($morning ? $date : $date + 1) :
        $date + 0.5;
    my $value = abs($try) > 1 ?
        sine_offset(dt_from_moment($alt), $location, $alpha) :
        $try;

    if (abs($value) <= 1) {
        return moment( local_from_apparent(
            dt_from_moment(
                $date + 0.5 +
                ($morning ? -1 : 1) *
                (mod(0.5 + asin_deg($value) / 360, 1) - 0.25)
            )
        ));
    } else {
        return undef;
    }
}

1;

__END__

=head1 NAME

DateTime::Util::Astro::Sun - Functions To Calculate Solnar Data

=head1 SYNOPSIS

  use DateTime::Util::Astro::Sun qw(solar_longitude);
  my $longitude = solar_longitude($dt);

=head1 DESCRIPTION

This module provides functions to calculate solar data, but its main
focus is to provide just enough functionality that allows us to
create lunisolar calendars and other DateTime related modules.

This module is a *straight* port from "Calendrical Calculations" [1] --
and therefore there are places where things can probably be "fixed" so
that they look more like Perl, as well as places where we could
leverage the DateTime functionalities better. If you see things that
doesn't quite look right (in Perl), that's probably because of that.

=head2 Notes On Accuracy

Before you use this module, please be aware that this module was originally
created B<solely> for the purpose of creating a lunisolar calendar for
the DateTime project (http://datetime.perl.org). 

We used [1] as basis for our calculations. While for most purposes the
results are accurate enough, you should note that the calculations from
this book are I<approximations>. 

Obviously we would like to make this module as good as possible, but
there's only so much you can do in the accuracy department. However, having
L<GMP|http://www.swox.com/gmp> and Math::BigInt::GMP may help a little bit.

This module by default uses Perl's arbitrary precision calculation module
Math::BigFloat. However, this adds a fair amount of overhead, and you will
see a noticeable difference in execution speed. This is true even if you
use GMP.

If you are willing to trade accuracy for speed, you may override the
class variable from DateTime::Util::Calc to toggle off the use of Math::BigFloat:

  use DateTime::Util::Astro::Sun qw(solar_longitude);
  local $DateTime::Util::Calc::NoBigFloat = 1;

  my $x = solnar_longitude($dt);

=head1 FUNCTIONS

=head2 solar_longitude($dt)

Given a DateTime object $dt, calculates the solar longitude at that time.

=head1 AUTHOR

Daisuke Maki E<lt>daisuke@cpan.orgE<gt>

=head1 REFERENCES

  [1] Edward M. Reingold, Nachum Dershowitz
      "Calendrical Calculations (Millenium Edition)", 2nd ed.
       Cambridge University Press, Cambridge, UK 2002

=head1 SEE ALSO

L<DateTime>
L<DateTime::Event::Lunar>
L<DateTime::Event::SolarTerm>
L<DateTime::Util::Astro::Common>
L<DateTime::Util::Astro::Moon>

=cut

