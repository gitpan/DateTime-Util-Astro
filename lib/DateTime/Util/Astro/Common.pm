# Please see file "LICENSE" for license information on code from
# "Calendrical Calculations".

package DateTime::Util::Astro::Common;
use strict;
use vars qw($VERSION @ISA @EXPORT_OK);

BEGIN
{
    $VERSION = '0.01';
    @ISA = qw(Exporter);
    @EXPORT_OK = qw(
        aberration
        dt_from_dynamical
        dynamical_moment_from_dt
        ephemeris_correction
        equation_of_time
        julian_centuries
        local_from_apparent 
        nutation
        obliquity
        standard_from_local
        standard_from_universal
        universal_from_local
        universal_from_standard
        SPRING SUMMER AUTUMN WINTER
        MEAN_TROPICAL_YEAR
        RD_MOMENT_1900_JAN_1
        RD_MOMENT_1810_JAN_1
        RD_MOMENT_J2000
    );
}

use DateTime;
use DateTime::Util::Calc
     qw(angle polynomial sin_deg cos_deg tan_deg bigfloat
        min moment dt_from_moment bf_downgrade);
use Math::Trig qw(pi);
use Params::Validate();

# I got the following from (DateTime->new(...)->utc_rd_values)[0]
use constant RD_MOMENT_1900_JAN_1   => 693596;
use constant RD_MOMENT_1810_JAN_1   => 660724.5;
use constant RD_MOMENT_J2000 => 730120.5;
use constant MEAN_TROPICAL_YEAR => 365.242189;
use constant SPRING => 0;
use constant SUMMER => 90;
use constant AUTUMN => 180;
use constant WINTER => 270;

# p.170
sub standard_from_universal
{
    my($dt, $location) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { isa => 'Astro::Earth::Location' }
    );
    dt_from_moment(moment($dt) + $location->zone / 24);
}

# p.170
sub universal_from_standard
{
    my($dt, $location) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { isa => 'Astro::Earth::Location' }
    );
    dt_from_moment(moment($dt) - $location->zone / 24);
}

# p.170
sub standard_from_local
{
    my($dt, $location) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { isa => 'Astro::Earth::Location' }
    );
    standard_from_universal(dt_from_moment(
        universal_from_local($dt, $location)), $location);
}

# p.170
sub local_from_standard
{
    my($dt, $location) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { isa => 'Astro::Earth::Location' }
    );
    local_from_universal(dt_from_moment(
        universal_from_standard($dt, $location)), $location);
}

# p.169
sub local_from_universal
{
    my($dt, $location) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { isa => 'Astro::Earth::Location' }
    );
    moment($dt) + $location->longitude / 360;
}
# p.169
sub universal_from_local
{
    my($dt, $location) = Params::Validate::validate_pos(@_,
        { isa => 'DateTime' }, { isa => 'Astro::Earth::Location' }
    );
    moment($dt) - $location->longitude / 360;
}

# [1] p172
sub dynamical_moment_from_dt
{
    my($dt) = Params::Validate::validate_pos(@_, { isa => 'DateTime' });
    return moment($dt) + ephemeris_correction($dt);
}

sub dt_from_dynamical
{
    my($t) = Params::Validate::validate_pos(@_, { type => Params::Validate::SCALAR() } );
    return dt_from_moment(
        $t - ephemeris_correction(dt_from_moment($t)));
}


# [1] p180
sub obliquity
{
    my $dt = shift;
    my $c  = julian_centuries($dt);
    return polynomial($c,
        angle(23, 26, 21.448),
        -1 * angle(0, 0, 46.8150),
        -1 * angle(0, 0, 0.00059),
        angle(0, 0, 0.001813)
    );
}

# [1] p171 + errata 158
sub ephemeris_correction
{
    my $dt = shift;

    # we need a gregorian calendar, so make sure $dt is just 'DateTime'
    if (ref($dt) ne 'DateTime') {
        $dt = DateTime->from_object(object => $dt);
    }

    my $year = $dt->year;

    # XXX - possible optimization

    my $c = EC_C($year);
    my $x = EC_X($year);

    my $correction;
    if (1988 <= $year && $year <= 2019) {
        $correction = EC1($year - 1933);
    } elsif (1900 <= $year && $year <= 1987) {
        $correction = EC2($c);
    } elsif (1800 <= $year && $year <= 1899) {
        $correction = EC3($c);
    } elsif (1700 <= $year && $year <= 1799) {
        $correction = EC4($year - 1700);
    } elsif (1620 <= $year && $year <= 1699) {
        $correction = EC5($year - 1600);
    } else {
        $correction = EC6($x);
    }

    return $correction;
}
sub EC_C {
    (bigfloat(RD_MOMENT_1900_JAN_1 -
        (DateTime->new(year => $_[0], month => 7, day => 1)->utc_rd_values)[0]))    / 36525;
}
sub EC_X {
    bigfloat(RD_MOMENT_1810_JAN_1 -
        (DateTime->new(year => $_[0], month => 1, day => 1)->utc_rd_values)[0]);
}
sub EC1 {
    $_[0] / (24 * 3600) }
sub EC2 {
    polynomial($_[0], -0.00002, 0.000297, 0.025184,
        -0.181133, 0.553040, -0.861938, 0.677066, -0.212591);
}
sub EC3 {
    polynomial($_[0],
        -0.00009, 0.003844, 0.083563, 0.865736, 4.867575,
        -15.845535, 31.332267, -38.291999, 28.316289,
        11.636204, 2.043794);
}
sub EC4 {
    polynomial($_[0], 8.118780842, -0.005092142,
        0.003336121, -0.0000266484) / (24 * 3600);
}
sub EC5 {
    polynomial($_[0], 196.58333,
        -4.0675, 0.0219167) / ( 24 * 3600 )
}
sub EC6 {
    (($_[0] ** 2 / 41048480 ) - 15) / ( 24 * 3600 ) }

# [1] p.183
sub aberration
{
    my $dt = shift;

    my $c = julian_centuries($dt);
    return 0.0000974 * cos_deg(177.63 + 35999.01848 * $c) - 0.0005575;
}

# [1] p.172
sub julian_centuries
{
    my($dt) = Params::Validate::validate_pos(@_, { isa => 'DateTime' });
    return (dynamical_moment_from_dt($dt) - RD_MOMENT_J2000) / 36525;
}

# [1] p.182
sub nutation
{
    my $dt = shift;

    my $c = julian_centuries($dt);
    my $A = polynomial($c, 124.90, -1934.134, 0.002063);
    my $B = polynomial($c, 201.11, 72001.5377, 0.00057);

    return bigfloat(-0.004778) * sin_deg($A) + 
        bigfloat(-0.0003667) * sin_deg($B);
}

# [1] p.177
sub equation_of_time
{
    my $dt = shift;
    my $c = julian_centuries($dt);
    my $longitude = polynomial($c,
        280.46645, 36000.76983, 0.0003032);
    my $anomaly = polynomial($c,
        357.52910, 35999.05030, -0.0001559, -0.00000048);
    my $eccentricity = polynomial($c,
        0.016708617, -0.000042037, -0.0000001236);
    my $epsilon = obliquity($dt);
    my $y = tan_deg($epsilon / 2) ** 2;

    my $equation = 
        (bigfloat(1) / (2 * pi)) *
        ($y * sin_deg(2 * $longitude) +
            (-2 * $eccentricity * sin_deg($anomaly)) + 
            (4 * $eccentricity * $y * sin_deg($anomaly) * cos_deg(2 * $longitude)) +
            (-0.5 * ($y ** 2) * sin_deg(4 * $longitude)) +
            (-1.25 * ($eccentricity ** 2) * sin_deg(2 * $anomaly)));

    my $sign = $equation >= 0 ? 1 : -1;
    $sign * min(abs($equation), 0.5);
}

# [1] p.178
sub local_from_apparent
{
    my $dt = shift;
    my $delta = equation_of_time($dt);
    if ($delta) {
        $dt->subtract(seconds => bf_downgrade($delta * 24 * 3600));
    }
    $dt;
}

#BEGIN
#{
#    if (eval { require Memoize } && !$@) {
#        Memoize::memoize('EC_C');
#        Memoize::memoize('EC_X');
#        Memoize::memoize('EC1');
#        Memoize::memoize('EC2');
#        Memoize::memoize('EC3');
#        Memoize::memoize('EC4');
#        Memoize::memoize('EC5');
#        Memoize::memoize('EC6');
#        Memoize::memoize('julian_centuries');
#    }
#}

package DateTime::Util::Astro::::Location;
use strict;
use Params::Validate ();

sub new
{
    my $class = shift;
    my %args = Params::Validate::validate(@_, {
        longitude => { type => Params::Validate::SCALAR(), default => 0 },
        latitude  => { type => Params::Validate::SCALAR(), default => 0 },
        elevation => { type => Params::Validate::SCALAR(), default => 0 },
        zone      => { type => Params::Validate::SCALAR(), default => 0 },
    });

    bless \%args, $class;
}

BEGIN
{
    foreach my $f qw(longitude latitude elevation zone) {
        eval sprintf( <<'EOM', $f, $f, $f);
sub %s
{
    my $self = shift;
    my $ret = $self->{%s};
    if (@_) {
        my $val = Params::Validate::validate_pos(@_,
            { type => Params::Validate::SCALAR() });
        $self->{%s} = $val;
    }
    return $ret;
}
EOM
        die if $@;
    }
}

1;

__END__

=head1 NAME

DateTime::Util::Astro::Common - Common Utilities For Astronomical Calendar Calculations

=head1 SYNOPSIS

  use DateTime::Util::Astro::Common qw(
    aberration
    dt_from_dynamical
    dynamical_moment_from_dt
    ephemeris_correction
    equation_of_time
    julian_centuries
    local_from_apparent 
    nutation
    obliquity
    standard_from_local
    standard_from_universal
    universal_from_local
    universal_from_standard
    SPRING
    SUMMER
    AUTUMN
    WINTER
    MEAN_TROPICAL_YEAR
    RD_MOMENT_1900_JAN_1
    RD_MOMENT_1810_JAN_1
    RD_MOMENT_J2000
  );

  my $location = DateTime::Util::Astro::Location->new(
    longitude => $longitude,
    latitude  => $latitude,
    zone      => $zone,
    elevation => $elevation
  );

=head1 DESCRIPTION

DateTime::Util::Astro::Location implements some functions that are commonly
used for astronomical calculations. As with other DateTime::Util::Astro::
modules this module only implements the bare minimum required to make
astronomical calendars.

=head1 FUNCTIONS

=head2 aberration($dt)

Calculates the effect of the sun's moving during the time its light takes
takes to reach the Earth

=head2 dt_from_dynamical($moment)

=head2 dynamical_moment_from_dt($dt)

=head2 ephemeris_correction($dt)

Calculates the offset from "dynamical time", which is caused by the
retarding effects of tide and other atmospheric conditions.

=head2 equation_of_time($dt)

Calculates the difference between "apparent midnight" and the "mean midnight"

=head2 julian_centuris($dt)

Calculates the fractional number of centuries since January 1, 2000 (Gregorian).

=head2 local_from_apparent($dt)

=head2 nutation($dt)

Calculates the effect caused by the wobble of the Earth.

=head2 obliquity($)

Calculates the inclination of th Earth

=head2 standard_from_local($dt,$location)

=head2 standard_from_universal($dt,$location)

=head2 universal_from_local($dt,$location)

=head2 universal_from_standard($dt,$location)

=head1 CONSTANTS

=head2 SPRING

=head2 SUMMER

=head2 AUTUMN

=head2 WINTER

The solar longitude for equinoxes and solstices.

=head2 MEAN_TROPICAL_YEAR

The time it takes for the mean sun to return to the same position relative
to the celestial equator

=head2 RD_MOMENT_1900_JAN_1

=head2 RD_MOMENT_1810_JAN_1

=head2 RD_MOMENT_J2000

=head1 TODO

DateTime::Util::Astro::Location probably isn't worth existing, as I believe
the math involved with locations can probably be done via DateTime itself,
or at least via DateTime::Locale and DateTime::TimeZone objects.

Hence it is foraseeable that DateTime::Util::Astro::Location will be
phased out eventually and I've opted to keep it as a private class within 
his module, so that it doesn't pollute users' file systems.

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
L<DateTime::Util::Astro::Moon>
L<DateTime::Util::Astro::Sun>

=cut


