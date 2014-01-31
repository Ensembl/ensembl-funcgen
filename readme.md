# ensembl-funcgen
The Ensembl Funcgen API (Application Programme Interface) serves as a middle layer between the underlying MySQL database and aims to encapsulate the database layout by providing high level access to the database. The team providing the data is Ensembl Regulation, but for historic reasons, the API and databases are called Funcgen.


## Usage

To clone and to bring onto your bash path use the following:

```
git clone https://github.com/Ensembl/ensembl-funcgen.git
```

```
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
    #-verbose => 1, #specificy verbose to see exactly what it loaded
);
```


```
```

## Contact us
http://www.ensembl.org/contact
Dev mailing list: http://lists.ensembl.org/mailman/listinfo/dev
