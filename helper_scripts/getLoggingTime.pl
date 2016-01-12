# subroutine for getting logging time
# this script can be used to add logging time to output of your scripts
# timestamp getLoggingTime()
sub getLoggingTime {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
        my $nice_timestamp = sprintf ( "%04d%02d%02d_%02d:%02d:%02d",
					$year+1900,$mon+1,$mday,$hour,$min,$sec);
	return $nice_timestamp;
}
1;
