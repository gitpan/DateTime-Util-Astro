#!perl

# data taken from http://aa.usno.navy.mil/data/docs/EarthSeasons.html

use strict;
use Test::More qw(no_plan);
BEGIN
{
    use_ok("DateTime::Util::Astro::Sun", qw(solar_longitude solar_longitude_after solar_longitude_before));
}

use constant ALLOWED_LONGITUDE_DELTA => 0.01;
use constant ALLOWED_DATE_SECOND_DELTA => 60 * 10;

my @dataset = (
    [  90, 1992,  6, 21,  3, 14 ],
    [ 270, 1992, 12, 21, 14, 43 ],
    [  90, 1993,  6, 21,  9, 00 ],
    [ 270, 1993, 12, 21, 20, 26 ],
    [  90, 1994,  6, 21, 14, 48 ],
    [  90, 1999,  6, 21, 19, 49 ]
);
foreach my $data (@dataset) {
    my $expected = $data->[0];
    my $dt    = DateTime->new(
        time_zone => 'UTC',
        year => $data->[1],
        month => $data->[2],
        day => $data->[3],
        hour => $data->[4],
        minute => $data->[5]
    );
    my $l     = solar_longitude($dt);
    my $delta = abs($expected - $l);
    ok($delta < ALLOWED_LONGITUDE_DELTA, "Expected $expected, got $l. Delta was $delta");

    my $dt2 = solar_longitude_after($dt->clone->subtract(days => 10), $expected);

    $delta = abs($dt2->epoch - $dt->epoch);
    ok($delta < ALLOWED_DATE_SECOND_DELTA, "Expected $dt, got $dt2. Delta was $delta");

    my $dt3 = solar_longitude_before($dt->clone->add(days => 10), $expected);
    $delta = abs($dt3->epoch - $dt->epoch);
    ok($delta < ALLOWED_DATE_SECOND_DELTA, "Expected $dt, got $dt3. Delta was $delta");
}