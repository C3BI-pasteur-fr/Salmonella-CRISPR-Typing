#!/usr/bin/perl 
# CRISPRNet Version 1.0
# Written by G.Guigon pour serveur www2.pasteur.fr
# use warnings;

use strict;
use CGI;
use IO::File;
use lib '/local/web/www.pasteur.fr/cgi-bin/genopole/PF8/crispr/sql' ;
use SQL::connection;


local (our ($database, $stylesheet )); 
local (our (%system));
local ( our  $q       = new CGI);
local (our ($page, $login, $pwd, $user, ));

#####################################################

&initiate;
&printpage;
$database->close();

#####################################################

#initialisation des variables
sub initiate {

	$page = $q->param('page');
	$page = 'index_page' if ( !$page );
	
	
	#variables "d'environement"
	$ENV{PATH} = '/bin:/usr/bin';
	
	if ( !$system{'host'} ) {
		$system{'host'} = 'fauve.sis.pasteur.fr';
	}
	if ( !$system{'port'} ) {
		$system{'port'} = 5435;
	}
	if ( !$system{'user'} ) {
		$system{'user'} = 'bbpe';
	}
	if ( !$system{'pwd'} ) {
		$system{'pwd'} = 'pasteur1';
	}
	if ( !$system{'dbname'} ) {
		$system{'dbname'} = 'salmcrispr';
	}
	if ( !$system{'scriptname'} ) {
		$system{'scriptname'} = 'crispr_2.pl';
	}
	if ( !$system{'description'} ) {
		$system{'description'} = 'Salmonella CRISPR Database';
	}
	
	#connection a la base et preparation des consultations
	$database = SQL::connection->new();
	$database->connect( $system{'host'}, $system{'port'}, $system{'dbname'}, $system{'user'}, $system{'pwd'} );

}

sub printpage {
	
	#feuille de style
	my $stylesheet = "/recherche/genopole/PF8/crispr/stylesheet.PF8db.css" ;
	my $footerfile = "/export/home/gguigon/mlva_web/footer.html";
	
	#definition du type de page (texte, html)
	print $q->header;
	print $q->start_html(
			-title    => "$system{'description'}",
			-encoding => 'UTF-8',
			-style    => { 'src' => $stylesheet }
		);
	&print_header;
 	print "<title>Institut Pasteur CRISPR Query Page</title>";
	if ( $page eq 'index_page' ) {
		&index_page;
	}

	elsif ($page eq 'check_log') {
		&check_log;
	}
	elsif ($page eq 'browse') {
		&browse;
	}
	elsif ($page eq 'search_seq' ) {
		&search_seq;
	}
	elsif ($page eq 'one_seq' ) {
		&one_seq;
	}
	elsif ($page eq 'one_seq_query' ) {
		&one_seq_query;
	}
	
	elsif ($page eq 'search_refDB' ) {
		&search_refDB;
	}
	print "<br><br>";
	&print_footer;
		
}

sub print_header {
print <<"HTML";
<body>
<table border="0" cellpadding="0" cellspacing="0" width="100%">
  <tbody>
    <tr>
      <td colspan="2" width="923"><!-- menu haut -->
        <table border="0" cellpadding="0" cellspacing="0" width="100%">
          <tbody>
            <tr>
              <td width="20%"><div align="center"> <img src="/recherche/genopole/PF8/crispr/images/logo%20pasteur.jpg" alt="" height="101" width="256" border="0"></div></td>
              <td width="60%"><div align="center"> <img src="/recherche/genopole/PF8/crispr/images/bandeau-accueil.jpg" alt="" height="89" width="481" border="0"></div></td>
              <td width="20%"><div align="right"> <img src="/recherche/genopole/PF8/crispr/images/pf8logo.jpg" alt="" height="122" width="141" border="0"></div></td>
            </tr>
          </tbody>
        </table></td>
    </tr>
	</table>
</body>
HTML
}

#Footer de toutes les pages/
sub print_footer {
print <<"HTML";
<body>
	<table border="0" cellpadding="0" cellspacing="0" width="100%">
	  <tbody>
	  	<tr>
	      <td colspan="3" align="center"  height="40">
	          <span class="cbblack"><b><a href="http://www.pasteur.fr/cgi-bin/genopole/PF8/crispr/crispr_2.pl"></b>Return to Home Page </a></span>
	      </td>
	    </tr>
	    <tr bgcolor="#364387">
	      <td width="917" height="10" colspan="3"></td>
	    </tr>
		<tr>
	      <td colspan="3" align="center"  height="40">
	          <i><span class="csblack">This site was developed at the technological platform <a href="http://www.pasteur.fr/recherche/genopole/PF8"><b>Genotyping of Pathogens and Public Health</b> </a>of Institut Pasteur by G.Guigon. </span></i>
	      </td>
	    </tr>
	  </tbody>
	</table>
</body>
HTML
}

#Page index pour authentification/
sub index_page {	
	print "<h1 class='Titre'>Welcome to the $system{'description'}</h1>\n\n";
	print "	<p>This server provides tools for querying the CRISPR database:</p>";		
	print "	<table class=\"formtable\" style=\"width:100%\"><tr><td>";
	print "	<p>To access to datasets, please log in:<br><br></p>";
 	print $q->start_form;
	print "<table><tr><td>Login</td><td>";
	print $q->textfield( -name => 'login', -size => '10' );
	print "</td></tr><tr><td>Password</td><td>";
	print $q->password_field( -name => 'password', -size => '10' );
	print "</td><td/><td>";
	$q->param( 'page', 'check_log' );
	print $q->hidden('page');
	print $q->hidden('form', 1);
	print $q->submit ( -name => 'Submit');
	print "</td></tr></table></td></tr></table>";
	print $q->end_form;
}

#Verification login et mot de passe avant d'accéder au menu/
sub check_log {
	if ($q->param('login')) {
		my $login = $q->param('login');
		my $pwd = $q->param('password');
		#if ($login eq $system{'user'} && $pwd eq $system{'pwd'}) {
		if ($login eq 'crisprteam' && $pwd eq 'febit$2008') {
			&menu_page;
		}
		else {
			&check_log;
			print "<br>Unknown user. Please enter correct login and password or go back to the <a href='http://www.pasteur.fr/cgi-bin/genopole/PF8/crispr/crispr_2.pl'>Home Page</a>";
		}
	}
}

#Page d'accueil, choix de la requete/
sub menu_page {	
print <<"HTML";
	<body>
	<h1 class='Titre'>Welcome to the $system{'description'}</h1>\n\n
	<br><br>
	<table><tr><td style="width: 50%; padding-right:3em; vertical-align: top">
	<p>This server provides tools for querying the CRISPR database:</p>	
	<br><br>
	<h2>Database queries:</h2>
	<br><br>
	<ul>
	<li><a href="$system{'scriptname'}?page=browse&amp">Browse spacers dictionary for CRISPR-1 and CRISPR-2 loci</a><br> - in 'spacer' experiment</li>
	<li><a href="$system{'scriptname'}?page=search_seq&amp">Browse spacers composition for all sequenced strains in database</a><br> - in 'CRISPR_concat1et2' experiment
	<li><a href="$system{'scriptname'}?page=search_refDB&amp">Browse spacers composition for the Reference strains</a><br>  - generate Binary Matrix of spacers for profile analysis.<br> - generate DR Distribution File
	<li><a href="$system{'scriptname'}?page=one_seq&amp">Search spacers composition for query</a><br>  
	- obtain a spacer profile for CRISPR-1 and CRISPR-2 loci from concatenate sequence in 'CRISPR_concat1et2' experiment.</li>
	</ul></table></body>
HTML
}#/

#Browse spacer composition for all the strains of database/
sub search_seq {
		my $out = IO::File->new();
		my $file = "/local/web/www.pasteur.fr/htdocs/recherche/genopole/PF8/tmp/Browse-spacers-composition.txt";
		chmod (0777, "/local/web/www.pasteur.fr/htdocs/recherche/genopole/PF8/tmp/Browse-spacers-composition.txt");
		$out->open(">$file") or die ("Impossible d'ouvrir le fichier : $!\n");
		print "<body>";
		print "<h1>Strain spacer composition for $system{'description'}: </h1>\n";
		print "<br><br>";
		print "<p align=center><a href=\"/recherche/genopole/PF8/tmp/Browse-spacers-composition.txt\"> Download the text-tabulated file</a></p>";
		print "<br>";
		my %serotype;
		my $td ='td1';
		#recuperation des serotypes pour chaque souche
		my @ser = $database->runqueries("SELECT DISTINCT \"ENTRYTABLE\".\"KEY\", \"ENTRYTABLE\".\"Serotype\"  FROM \"ENTRYTABLE\"");
		foreach my $ser (@ser) {
			my @sero = split('///',$ser); 
			keys %serotype = $sero[0];
		    $serotype{$sero[0]} = $sero[1];
		    }
		print "<table class='resultstable'>";
		#recuperation des sequences pour chaque souche et identification des spacers
		my @strains = $database->runqueries("SELECT DISTINCT \"SEQUENCES\".\"KEY\", \"SEQUENCES\".\"SEQUENCE\"  FROM \"SEQUENCES\" WHERE \"SEQUENCES\".\"EXPERIMENT\"=\'CRISPR_concat1et2\'");
		foreach my $strains (@strains) {
			my @seq = split('///',$strains);
			#my $seq = 'TGCTGTTGAAACGTGTTTATCCCCGCTGGCGCGGGGAACACACTTCGAGTGATGGAAAGACAGGAGTTGTCTTCGGTTTATCCCCGCTGGCGCGGGGAACACACTTTACATCATTAGCTATGTCCAAATCATACCGGTTTATCCCTGCTGGCGCGGGGAACACAGCATTGGCGTTTTCCGTTGCGTCAATGGCTTCGGTTTATCCCCGCTGACGCGGGGAACACACCCCTTTAATCTGTTCATTGCTCCAGCCACGCGGTTTATCCTCGCTGGCGCGGGGAACACAAATTGCAATCCCCCCACGCCACTGGTCAGACCGGTTTATCCCCGCTGGCGCGGGGAACACGCAGGATGACCGGCGGCTGGTTCTCCTGCGTTCGGTTTATCCCCGCTGGCGCGGGGAACACGATGAGCAACACGCCCGCACTGGCGTAACTTACGGTTTATCCCCGCTGGCGCGGGGAACACATATCCAGCCTGCGATAGAGATAGTTTTCCGGCGGTTTATCCCCGCTGGCGCGGGGAACACTTTTCAGCTCGCCGCAGAATATCAGTGAGGCGCGGTTTATCCCCGCTGGCGCGGGGAACACTTTCTCGACCTCGTCGCCGTAACGTTCGGCCGCGGTTTATCCCCGCTGGCGCGGGGAACACCAGGTGCTCGCTATGATGTTTATCTTGATGGGCGGTTTATCCCCGCTGGCGCGGGGAACACTGCTGCCGTGGAAGAGTTTCATTTGATAAATCCGGTTTATCCCCGCTGGCGCGGGGAAC|AAGCACTGGACATTATTTTATCGCGTTGTGGCAACCGGACAGTGCATTTGTTGGATCCGCAACACGGTTGATGATGCGCTGGATACTTACCAGCAGTTGCTACACGAGGGAATCGTTCCACAGCAGGATCTTTTGCTTTTCCACAGCCGTTTTGCTTTTATCGATCGTATTGCTATTGAAAATAAAACGTTAAACTGGTTTGGTAACAATGCCCCCGTTTCTGAACGACGTGGTAAGGTATTAATCGCCACACAAGTCGTTGAACAAAGTCTCGATCTGGACTTTGACTGGATGATCACCGATCTGGCGCCGATTGATTTATTGATCCAGCGCGCCGGTCGCCTGCAACGCCATATCCGCGATGCCCATGGGCAACGAAAATCTACGCTTCCGGATGAACGCCAGCCACCAATATTGCATATTCTGGCCCCTCACTGGCAGGAGCAGGCAGAAGAGAGCTGGTTGGGGCAAGAGCTAAAAGGTACCGGCTTTGTTTACACCGATCATGCCTGCCTGTGGCGTACACAGGCGCTGCTGCGTCAACATGGTGAGATCCGTATGCCGGATAATGCCAGAGCACTCGTTGACGGCGTATACGAGCAGAAAATTGCTGCGCCAGCAGGCTTGCAGACGATCTCTGACGTTGCCTTTGGCAAAGTCCTTAGCCAACGTTCAGTGGCTGCGCAAAATTTATTACGTTAC|TGCCAATCGGCAGGGGGTCACGCGTGTACATCATTCGGTAGCGTGAATAGTCTACCTCACACGACAATTTTTTTATTTTTCGGGTAATACCGAGAAAAGTCATTGCCGCATAACGTTCGCGCGCAATATCGGGATGGCAGATGATTCGACTGTTATCAGGAAGCCACGGCACGCCGCCGCAA';
			my $new = $seq[1];
			if (defined $new){
				my @spacer = $database->runqueries("SELECT DISTINCT \"SEQUENCES\".\"KEY\", \"SEQUENCES\".\"SEQUENCE\"  FROM \"SEQUENCES\" WHERE \"SEQUENCES\".\"EXPERIMENT\"=\'spacer\'");
				foreach my $spac (@spacer) {
					my @tab = split('///',$spac);
		            #print "<table class='resultstable'><tr class='td1'><td>$tab[0]</td><td>$spacers{$tab[0]}</td></tr>";
					if ($new =~/^.*$tab[1].*$/i) {
						$new =~s/$tab[1]/-$tab[0]/gi;#/
					}
		        }
			    if ($td eq 'td1') {
			    	$td = 'td2';
			    	}
			    else {$td = 'td1';}
				print "<tr class=$td><td><b>$seq[0]: </b><td><b>$serotype{$seq[0]}</b></td><td align=left>$new</td></tr>";
				print $out $seq[0].":"."\t".$serotype{$seq[0]}."\t".$new."\n";
			}
		}
		print "</table></body>";
		$out ->close();
}
#/
#Browse spacer composition for a query, paste page/
sub one_seq {
		print "<body>";
		print "<h1>Strain spacer composition for $system{'description'}: </h1>\n";
		print "<br><br>";
		print $q->start_form;
		$q->param( 'page', 'one_seq_query' );
		print $q->hidden('page');
		print "<p>Enter sequence (DNA) below using copy and paste:</p>\n";
		print $q->textarea(
			-name     => 'sequence',
			-rows     => 8,
			-columns  => 60,
			-override => 1
		);
		print "<p>\n";
		#print $q->hidden('form', 1);
		print $q->reset;
		print $q->submit;
		print "</p>\n";
		print $q->endform;
		print "</body>";
}

#Browse spacer composition for the reference strains of database and create the binary file/
sub search_refDB {
		my $out = IO::File->new();
		my $file = "/local/web/www.pasteur.fr/htdocs/recherche/genopole/PF8/tmp/Binary-matrix.txt";
		chmod (0777, "/local/web/www.pasteur.fr/htdocs/recherche/genopole/PF8/tmp/Binary-matrix.txt");
		$out->open(">$file") or die ("Impossible d'ouvrir le fichier : $!\n");
		my $DR = IO::File->new();
		my $file2 = "/local/web/www.pasteur.fr/htdocs/recherche/genopole/PF8/tmp/DR-distribution.txt";
		chmod (0777, "/local/web/www.pasteur.fr/htdocs/recherche/genopole/PF8/tmp/DR-distribution.txt");
		$DR->open(">$file2") or die ("Impossible d'ouvrir le fichier : $!\n");
		print "<body>";
		print "<h1>Strain spacer composition for $system{'description'}: </h1>\n";
		print "<br><br>";
		print "<p align=center><a href=\"/recherche/genopole/PF8/tmp/Binary-matrix.txt\"> Download the Binary Matrix File</a></p>";
		print "<p align=center><a href=\"/recherche/genopole/PF8/tmp/DR-distribution.txt\"> Download the DR distribution File</a></p>";
		print "<br>";
		my %serotype;
		my $td ='td1';
		print "<table class='resultstable'>";
		
		#recuperation des serotypes pour chaque souche
		my @ser = $database->runqueries("SELECT DISTINCT \"ENTRYTABLE\".\"KEY\", \"ENTRYTABLE\".\"Serotype\"  FROM \"ENTRYTABLE\"");
		foreach my $ser (@ser) {
			my @sero = split('///',$ser); 
			keys %serotype = $sero[0];
		    $serotype{$sero[0]} = $sero[1];
		    }
		
		#recuperation des sequences des spacers dans l'experiment refDB
		my @spacer = $database->runqueries("SELECT DISTINCT \"SEQUENCES\".\"KEY\", \"SEQUENCES\".\"SEQUENCE\"  FROM \"SEQUENCES\" WHERE \"SEQUENCES\".\"EXPERIMENT\"=\'spacer\'");

		#recuperation des sequences pour chaque souche et identification des spacers
		my @strains = $database->runqueries("SELECT DISTINCT \"SEQUENCES\".\"KEY\", \"SEQUENCES\".\"SEQUENCE\"  FROM \"SEQUENCES\" WHERE \"SEQUENCES\".\"EXPERIMENT\"=\'refDB\'");
		foreach my $strains (@strains) {
			my @seq = split('///',$strains);
			my $new = $seq[1];
			if (defined $new){
				foreach my $spac (@spacer) {
					my @tab = split('///',$spac);
		            #print "<table class='resultstable'><tr class='td1'><td>$tab[0]</td><td>$spacers{$tab[0]}</td></tr>";
					if ($new =~/^.*$tab[1].*$/i) {
						$new =~s/$tab[1]/-$tab[0]/gi;#/
					}
		        }
			    if ($td eq 'td1') {
			    	$td = 'td2';
			    	}
			    else {$td = 'td1';}
				print "<tr class=$td><td><b>$seq[0]: </b><td><b>$serotype{$seq[0]}</b></td><td align=left>$new</td></tr>";
				
				
			}
		}
		print "</table></body>";
		# creation du fichier binaire et de la distribution des DR
		print $out " "."\t";
		print $DR " "."\t";
		foreach my $strains (@strains) {
				my @seq1 = split('///',$strains);
				print $out $seq1[0]."(".$serotype{$seq1[0]}.")\t";
				print $DR $seq1[0]."(".$serotype{$seq1[0]}.")\t";
				}
		print $out "\n";
		print $DR "\n";
		foreach my $spac (@spacer) {
			my @tab = split('///',$spac);
			if ($tab[0]=~/DR.*/i){
				print $DR $tab[0]."\t";
				foreach my $strains (@strains) {
					my @seq2 = split('///',$strains);
					if ($seq2[1] =~/^.*$tab[1].*$/i){
					print $DR "1\t";
					}
					else {
					print $DR "0\t";
					}
				}
			print $DR "\n";
			}
			else {
			print $out $tab[0]."\t";
				foreach my $strains (@strains) {
					my @seq3 = split('///',$strains);
					if ($seq3[1] =~/^.*$tab[1].*$/i){
					print $out "1\t";
					}
					else {
					print $out "0\t";
					}
				}
			print $out "\n";
			}
		}		
		
		$out ->close();
}

#Browse spacer composition for a query, result page/
sub one_seq_query {
		print "<body>";
		print "<h1>Strain spacer composition for $system{'description'}: </h1>\n";
		print "<br><br>";
		if ($q->param('sequence')) {
		my $seq = $q->param('sequence');
			if (defined $seq){
				my @spacer = $database->runqueries("SELECT DISTINCT \"SEQUENCES\".\"KEY\", \"SEQUENCES\".\"SEQUENCE\"  FROM \"SEQUENCES\" WHERE \"SEQUENCES\".\"EXPERIMENT\"=\'spacer\'");
				foreach my $spac (@spacer) {
					my @tab = split ('///', $spac);
					if ($seq =~/^.*$tab[1].*$/i) {
						$seq =~s/$tab[1]/-$tab[0]/gi;#/
						}
				}
			print "<table class='resultstable'><tr class='td1'><td><b>Query: </b></td><td>$seq</td></tr></table></body>";
			}
		}
}
#"
sub composition () {
	my ($seq) = @_;
	my %occur;
	my $num_GC;
		foreach my $base (split('',$seq)) {
		$occur{$base}++;
		if(($base=~"G" || $base=~"g") or ($base=~"C" || $base=~"c")) #if it matches G or C increase counter
			{$num_GC++;
			}
		}
	my $GC_content=($num_GC/length($seq))*100; #/
	$occur{'GC'}=$GC_content;
	return (%occur);
}

#Browse spacer dictionary for all database/
sub browse {
		print "<body>";
		print "<h1>Browse $system{'description'}: </h1>\n";
		print "<br><br>";
		#my $sql = qq|"SELECT DISTINCT "SEQUENCES"."SEQUENCE" FROM "SEQUENCES" WHERE "SEQUENCES"."EXPERIMENT"='CRISPR-1'"|;
		#my @sets = $database->runqueries($sql);
		print "<h2>Distinct sequences of spacers for CRISPR-1 and CRISPR-2 loci</h2>";
		print "<table class='resultstable'><tr align=left><th>ID</th><th>Key</th><th>GC Content</th><th>Sequence</th></tr>";
		my @sets = $database->runqueries("SELECT DISTINCT \"SEQUENCES\".\"KEY\", \"SEQUENCES\".\"SEQUENCE\"  FROM \"SEQUENCES\" WHERE \"SEQUENCES\".\"EXPERIMENT\"=\'spacer\'");
				my $i;
				my $sum=0;
				foreach my $sets (@sets) {
					my @tab = split ('///', $sets);
					my (%occ) = &composition($tab[1]);
					$i++;
					$sum = $sum + $occ{'GC'};
					print "<tr align=left><td>$i</td><td>$tab[0]</td><td>$occ{'GC'}</td><td>$tab[1]</td></tr>";
				}
				my $mean=$sum/$i;
				print "</table>";
				print "<p>GC content average : $mean</p>";
				print "</body>";
}