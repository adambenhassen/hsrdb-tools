CREATE TABLE DeprecatedPatterns(ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, OldReferenceCode TEXT UNIQUE NOT NULL, NewReferenceCode TEXT NOT NULL );
CREATE TABLE literature_map(doc INTEGER PRIMARY KEY NOT NULL, pid INTEGER NOT NULL );
CREATE TABLE PatternPhase(ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, ReferenceCode TEXT NOT NULL, PhaseName TEXT NOT NULL, UNIQUE(ReferenceCode, PhaseName) );
CREATE TABLE Phases(ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, Name TEXT NOT NULL UNIQUE, Phase BLOB NOT NULL, Peaks BLOB );
CREATE TABLE PatternScan(ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, ReferenceCode TEXT NOT NULL, ScanName TEXT NOT NULL, UNIQUE(ReferenceCode, ScanName) );
CREATE TABLE Scans(ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, Name TEXT NOT NULL UNIQUE, Scan BLOB NOT NULL );
CREATE TABLE Properties(ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, Key TEXT NOT NULL UNIQUE, Value TEXT );
CREATE TABLE SubFileRef(pid INTEGER NOT NULL, SubFileID INTEGER NOT NULL );
CREATE TABLE Subfiles(ID INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, SubFile TEXT UNIQUE NOT NULL );
CREATE TABLE general_indexed(id INTEGER PRIMARY KEY NOT NULL, ProductID TEXT UNIQUE NOT NULL, XTSLSYS TEXT NOT NULL DEFAULT '', STATUS INTEGER, A REAL, B REAL, C REAL, ALPHA REAL, BETA REAL, GAMMA REAL, A_OVER_C REAL, Z INTEGER, DM REAL, DX REAL, D REAL, QUALFINAL TEXT NOT NULL DEFAULT '', IIC REAL, NumberOfElements INTEGER, XTLASPECTN INTEGER, XTLVOL REAL, CTIME INTEGER, MTIME INTEGER);
CREATE TABLE general_stored(pid INTEGER PRIMARY KEY NOT NULL, CELLED TEXT, SPGR TEXT, SPGRED TEXT, ANX REAL, XTLSG TEXT, XTLSGED TEXT, Comment TEXT, CAS TEXT, Lines BLOB NOT NULL, HKL BLOB, LinesI BLOB NOT NULL );
CREATE TABLE general_used(  pid INTEGER PRIMARY KEY NOT NULL,   ProductID TEXT UNIQUE NOT NULL );
CREATE TABLE pattern_elements (   pid INTEGER PRIMARY KEY NOT NULL,   e1 INTEGER,   e2 INTEGER );
CREATE TABLE pattern_strongestlines(  pid INTEGER PRIMARY KEY NOT NULL,   ProductID TEXT UNIQUE NOT NULL,   Lines BLOB NOT NULL );
CREATE INDEX literature_id ON literature_map(pid);
CREATE INDEX General_a ON general_indexed(a);
CREATE INDEX General_a_over_c ON general_indexed(a_over_c);
CREATE INDEX General_alpha ON general_indexed(alpha);
CREATE INDEX General_b ON general_indexed(b);
CREATE INDEX General_beta ON general_indexed(beta);
CREATE INDEX General_c ON general_indexed(c);
CREATE INDEX General_d ON general_indexed(d);
CREATE INDEX General_gamma ON general_indexed(gamma);
CREATE INDEX General_iic ON general_indexed(iic);
CREATE INDEX General_numberofelements ON general_indexed(numberofelements);
CREATE INDEX General_xtlaspectn ON general_indexed(xtlaspectn);
CREATE INDEX General_xtlvol ON general_indexed(xtlvol);
CREATE INDEX General_xtslsys ON general_indexed(xtslsys);
CREATE INDEX General_z ON general_indexed(z);
CREATE INDEX SubFileRef_ID ON SubFileRef(pid,subfileid);
CREATE INDEX SubFileRef_SubFileID ON SubFileRef(SubFileID, pid);
CREATE VIRTUAL TABLE literature_text using fts3(volume TEXT NOT NULL DEFAULT '', pages TEXT NOT NULL DEFAULT '', year TEXT NOT NULL DEFAULT '', authors TEXT NOT NULL DEFAULT '', journal TEXT NOT NULL DEFAULT '', referencetype Integer NOT NULL DEFAULT 0 )
/* literature_text(volume,pages,year,authors,journal,referencetype) */;
CREATE VIRTUAL TABLE general_text using fts3(  chemicalformula TEXT,   compoundname TEXT,   mineralname TEXT,   commonname TEXT,   empiricalformula TEXT DEFAULT '',   color TEXT )
/* general_text(chemicalformula,compoundname,mineralname,commonname,empiricalformula,color) */;
CREATE VIEW general as select id, ProductID, XTSLSYS, STATUS, A, B, C, ALPHA, BETA, GAMMA, A_OVER_C, Z, DM, DX, D, QUALFINAL, IIC, NumberOfElements, XTLASPECTN, XTLVOL, CTIME, MTIME, CELLED, SPGR, SPGRED, ANX, XTLSG, XTLSGED, Comment, CAS, Lines, HKL, LinesI, chemicalformula, compoundname, mineralname, commonname, empiricalformula, color, e1, e2 from general_indexed gi, general_stored gs, general_text gt, pattern_elements pe where id = gs.pid and id = gt.docid and id = id and id = pe.pid
/* general(id,ProductID,XTSLSYS,STATUS,A,B,C,ALPHA,BETA,GAMMA,A_OVER_C,Z,DM,DX,D,QUALFINAL,IIC,NumberOfElements,XTLASPECTN,XTLVOL,CTIME,MTIME,CELLED,SPGR,SPGRED,ANX,XTLSG,XTLSGED,Comment,CAS,Lines,HKL,LinesI,chemicalformula,compoundname,mineralname,commonname,empiricalformula,color,e1,e2) */;
CREATE TRIGGER general_insert INSTEAD OF INSERT ON general FOR EACH ROW BEGIN INSERT INTO general_indexed (ID, ProductID, XTSLSYS, STATUS, A, B, C, ALPHA, BETA, GAMMA, A_OVER_C, Z, DM, DX, D, QUALFINAL, IIC, NumberOfElements, XTLASPECTN, XTLVOL, CTIME, MTIME ) VALUES (new.id, new.ProductID, new.XTSLSYS, new.STATUS, new.A, new.B, new.C, new.ALPHA, new.BETA, new.GAMMA, new.A / new.C, new.Z, new.DM, new.DX, coalesce(nullif(max(new.DM, 0), 0), new.DX), new.QUALFINAL, new.IIC, new.NumberOfElements, new.XTLASPECTN, new.XTLVOL, new.CTIME, new.MTIME); INSERT INTO general_stored (pid, CELLED, SPGR, SPGRED, ANX, XTLSG, XTLSGED, Comment, CAS, Lines, HKL, LinesI) VALUES (new.id, new.CELLED, new.SPGR, new.SPGRED, new.ANX, new.XTLSG, new.XTLSGED, new.Comment, new.CAS, new.Lines, new.HKL, new.LinesI); INSERT INTO general_text (docid, chemicalformula, compoundname, mineralname, commonname, empiricalformula, color) VALUES (new.id, new.chemicalformula, new.compoundname, new.mineralname, new.commonname, new.empiricalformula, new.color); INSERT OR IGNORE INTO general_used (pid, ProductID) VALUES (new.id, new.ProductID); INSERT INTO pattern_elements (pid, e1, e2) VALUES (new.id, new.e1, new.e2); INSERT INTO pattern_strongestlines (pid, ProductID, Lines) VALUES (new.id, new.ProductID, substr(new.Lines, 1, 60)); UPDATE properties SET Value = (SELECT Value FROM properties WHERE Key = 'NoOfPatterns') + 1 WHERE Key='NoOfPatterns'; END;
CREATE TRIGGER general_update INSTEAD OF UPDATE ON general FOR EACH ROW BEGIN UPDATE general_indexed SET ProductID = new.ProductID, XTSLSYS = new.XTSLSYS, STATUS = new.STATUS, A = new.A, B = new.B, C = new.C, ALPHA = new.ALPHA, BETA = new.BETA, GAMMA = new.GAMMA, A_OVER_C = new.A / new.C, Z = new.C, DM = new.DM, DX = new.DX, D = coalesce(nullif(max(new.DM, 0), 0), new.DX), QUALFINAL = new.QUALFINAL, IIC = new.IIC, NumberOfElements = new.NumberOfElements, XTLASPECTN = new.XTLASPECTN, XTLVOL = new.XTLVOL, CTIME = new.CTIME, MTIME = new.MTIME WHERE ID = old.id; UPDATE general_stored SET CELLED = new.CELLED, SPGR = new.SPGR, SPGRED = new.SPGRED, ANX = new.ANX, XTLSG = new.XTLSG, XTLSGED = new.XTLSGED, Comment = new.Comment, CAS = new.CAS, Lines = new.Lines, HKL = new.HKL, LinesI = new.LinesI WHERE pid = old.id; UPDATE general_text SET chemicalformula = new.chemicalformula, compoundname = new.compoundname, mineralname = new.mineralname, commonname = new.commonname, empiricalformula = new.empiricalformula, color = new.color WHERE docid = old.id; UPDATE pattern_elements SET e1 = new.e1, e2 = new.e2 WHERE pid=old.id; UPDATE pattern_strongestlines SET ProductID = new.ProductID, Lines = substr(new.Lines, 1, 60) WHERE pid = old.id; END;
CREATE TRIGGER general_delete INSTEAD OF DELETE ON general FOR EACH ROW BEGIN DELETE FROM general_indexed WHERE ID = old.id; DELETE FROM general_stored WHERE pid = old.id; DELETE FROM general_text WHERE docid = old.id; DELETE FROM general_stored WHERE pid = old.id; DELETE FROM pattern_elements WHERE pid = old.id; DELETE FROM pattern_strongestlines WHERE pid = old.id; DELETE FROM literature WHERE pid = old.id; DELETE FROM SubfileRef WHERE pid = old.id; UPDATE properties SET Value = (SELECT Value FROM properties WHERE key = 'NoOfPatterns') - 1 WHERE Key = 'NoOfPatterns'; END;
CREATE VIEW literature AS SELECT doc, pid, volume, pages, year, authors, journal, referencetype FROM literature_map, literature_text WHERE literature_map.doc = literature_text.docid
/* literature(doc,pid,volume,pages,year,authors,journal,referencetype) */;
CREATE TRIGGER literature_insert INSTEAD OF INSERT ON literature FOR EACH ROW BEGIN INSERT INTO literature_map(doc, pid) VALUES (new.doc, new.pid); INSERT INTO literature_text(docid, volume, pages, year, authors, journal, referencetype) VALUES (new.doc, new.volume, new.pages, new.year, new.authors, new.journal, new.referencetype); END;
CREATE TRIGGER literature_delete INSTEAD OF DELETE ON literature FOR EACH ROW BEGIN DELETE FROM literature_text WHERE docid IN (SELECT doc FROM literature_map WHERE pid = old.pid); DELETE FROM literature_map WHERE pid = old.pid; END;
