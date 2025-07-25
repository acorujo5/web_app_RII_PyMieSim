from collections import namedtuple, OrderedDict
import os
import yaml
import sqlite3
import numpy as np

from refractiveindexmaster.refractiveindex2.refractivesqlite import material
from refractiveindexmaster.refractiveindex2.refractivesqlite.material import Material

Shelf = namedtuple('Shelf', ['shelf', 'name'])
Book = namedtuple('Book', ['book', 'name'])
Page = namedtuple('Page', ['page', 'name', 'path'])
Entry = namedtuple('Entry', ['id', 'shelf', 'book', 'page'])

_riiurl = "https://refractiveindex.info/download/database/" +\
          "rii-database-2019-02-11.zip"


class Database:
    def __init__(self, sqlitedbpath):
        '''
        Construct a database instance

        :param sqlitedbpath: The path of the sqlitedatabse
                             it has to exist even if you want to create
                             a new database
        '''
        self.db_path = sqlitedbpath
        if not os.path.isfile(sqlitedbpath):
            print("Database file not found.")
        else:
            print("Database file found at", sqlitedbpath)

    def create_database_from_folder(self, yml_database_path,
                                    interpolation_points=100):
        '''
        Create a sql database from a yml database path

        :param yml_database_path: The path to the yaml database
        :param interpolation_points: The number of interpolation_points to use
        '''
        create_sqlite_database(yml_database_path,
                               self.db_path,
                               interpolation_points=interpolation_points)

    def create_database_from_url(self,
                                 riiurl=_riiurl,
                                 interpolation_points=100):
        '''
        Create a sqlite database from an url

        :param riiurl: The url where to download the zip compressed
                       refractive index database from
        :param interpolation_points: The number of interpolation_points to use
        '''
        Database.DownloadRIIzip(riiurl=riiurl)
        self.create_database_from_folder(
            "database", interpolation_points=interpolation_points)

    def check_url_version(self):
        print(_riiurl)

    def get_all_shelves(self):

        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        query = "SELECT DISTINCT shelf FROM pages;"
        c.execute(query)

        results = [row[0] for row in c.fetchall()]
        conn.close()
        return results

    def get_books_in_shelf(self, shelf):

        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        query = f'SELECT DISTINCT book FROM pages WHERE shelf = "{shelf}"';
        c.execute(query)

        results = [row[0] for row in c.fetchall()]
        conn.close()
        return results

    def get_author(self, shelf, book):

        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        query = f'SELECT * FROM pages WHERE shelf = "{shelf}" AND book = "{book}"';
        c.execute(query)

        results = [row for row in c.fetchall()]

        names = [row[3] for row in results]
        conn.close()
        return names

    def search_custom(self, sqlquery):
        '''
        Make a custom sql query

        :param sqlquery: The sql query to make
        :retrurn: Return all results of the query
        '''
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute(sqlquery)
        results = c.fetchall()
        if len(results) == 0:
            print("No results found.")
        else:
            print(len(results), "results found.")
        conn.close()
        return results

    def search_pages(self, term="", exact=False):
        '''
        Search for pages by a looking for the searchterm
        in shelf, book, page and filename

        :param term: The search term to look for
        :param exact: If false term is extended to %term%
        '''
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        if not exact:
            c.execute('SELECT * FROM pages WHERE shelf like ? or book like'
                      '? or page like ? or filepath like ?',
                      ["%"+term+"%" for i in range(4)])
        else:
            c.execute('SELECT * FROM pages WHERE shelf like ? or book like'
                      '? or page like ? or filepath like ?',
                      [term for i in range(4)])
        results = c.fetchall()
        if len(results) == 0:
            print("No results found.")
        else:
            pass
            #print(len(results), "results found.")
            #columns = self._get_pages_columns()
            #print("\t".join(columns))
            #for r in results:
            #    print("\t".join(map(str, r[:])))
        conn.close()
        return results

    def search_id(self, pageid):
        '''
        Print page informations

        :param pageid: The id of the page to print
        '''
        info = self._get_page_info(pageid)
        if info is None:
            print("PageID not found.")
        else:
            print("\t".join(info.keys()))
            print("\t".join(map(str, info.values())))

    def search_n(self, n, delta_n):
        '''
        Search for materials with a fraction index between
        n and delta_n

        :param n: The lower bound of the fraction index
        :param delta_m: The upper bound of the fraction index
        '''
        print("*Search n =", n, "delta_n = ", delta_n)
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        interval = [n-delta_n, n+delta_n]
        c.execute('''select r.pageid,shelf,book,page,r.wave,r.refindex
                    from refractiveindex r join pages p on r.pageid = p.pageid
                    where refindex between ? and ? ''', interval)
        results = c.fetchall()
        if len(results) == 0:
            print("No results found.")
        else:
            print(len(results), "results found.")
            print("pageid|shelf|book|page|wavelength|n")
            for r in results:
                print(r)
        conn.close()

    def search_k(self, k, delta_k):
        '''
        Search for materials with an extinction coefficient between
        k and delta_k

        :param k: The lower bound of the extinction coefficient
        :param delta_k: The upper bound of the extinction coefficient
        '''
        print("*Search k =", k, "delta_k = ", delta_k)
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        interval = [k-delta_k, k+delta_k]
        c.execute('''select e.pageid,shelf,book,page,e.wave,e.coeff
                    from extcoeff e join pages p on e.pageid = p.pageid
                    where coeff between ? and ?''', interval)
        results = c.fetchall()
        if len(results) == 0:
            print("No results found.")
        else:
            print(len(results), "results found.")
            print("pageid|shelf|book|page|wavelength|k")
            for r in results:
                print(r)
        conn.close()

    def search_nk(self, n, delta_n, k, delta_k):
        '''
        Search for materials with fraction indice and extinction
        coefficient between n and delta_n and k and delta_k

        :param n: The lower bound of the fraction index
        :param delta_n: The upper bound of the fraction index
        :param k: The lower bound of the extinction coefficient
        :param delta_k: The upper bound of the extinction coefficient
        '''
        print("*Search n =", n, "delta_n = ", delta_n, "k = ", k,
              "delta_k = ", delta_k)
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        interval = [n-delta_n, n+delta_n, k-delta_k, k+delta_k]
        c.execute('''select r.pageid, shelf, book, page, r.wave, r.refindex,
                     e.coeff from refractiveindex r join extcoeff e on
                     r.pageid = e.pageid and r.wave = e.wave join pages p
                     on r.pageid = p.pageid where refindex between
                     ? and ? and coeff between ? and ?''', interval)
        results = c.fetchall()
        if len(results) == 0:
            print("No results found.")
        else:
            print(len(results), "results found.")
            print("pageid|shelf|book|page|wavelength|n|k")
            for r in results:
                print(r)
        conn.close()

    def get_material(self, pageid):
        '''
        Get the material from a pageid

        :param pageid: The pageid of the material
        :returns: Material
        '''
        pagedata = self._get_page_info(pageid)
        if pagedata is None:
            print("PageID not found.")
            return None
        else:
            conn = sqlite3.connect(self.db_path)
            c = conn.cursor()
            wavelengths_r = None
            wavelengths_e = None
            refractive = None
            extinction = None
            if pagedata['hasrefractive'] == 1:
                c.execute('''select wave,refindex
                            from refractiveindex
                            where pageid = ?
                            order by wave asc''', [pageid])
                results = c.fetchall()
                wavelengths_r = [r[0] for r in results]
                refractive = [r[1] for r in results]
            if pagedata['hasextinction'] == 1:
                c.execute('''select wave,coeff
                            from extcoeff
                            where pageid = ?
                            order by wave asc''', [pageid])
                results = c.fetchall()
                wavelengths_e = [r[0] for r in results]
                extinction = [r[1] for r in results]
            conn.close()
            print("Material", pagedata['filepath'], "loaded.")
            return Material.FromLists(pagedata,
                                      wavelengths_r=wavelengths_r,
                                      refractive=refractive,
                                      wavelengths_e=wavelengths_e,
                                      extinction=extinction)

    def get_material_n_numpy(self, pageid):
        '''
        Get the refraction index of a material

        :param pageid: The pageid of the material
        :return: The refraction data as a numpy array
        '''
        mat = self.get_material(pageid)
        if mat is None:
            return None
        n = mat.get_complete_refractive()
        if n is None:
            print("Material has no refractive data.")
            return None
        return np.array(n)

    def get_material_k_numpy(self, pageid):
        '''
        Get the extinction coefficient of a material

        :param pageid: The pageid of the material
        :return: The extinction data as a numpy array
        '''
        mat = self.get_material(pageid)
        if mat is None:
            return None
        k = mat.get_complete_extinction()
        if k is None:
            print("Material has no extinction data.")
            return None
        return np.array(k)

    def get_material_csv(self, pageid, output="", folder=""):
        '''
        Safe a material as a comma seperated value list

        :param pageid: The pageid of the material
        :param output: The name of the output file
                       default is [pageid][shelf][book][page].csv
        :param folder: The output folder default is ./
        '''
        mat = self.get_material(pageid)
        if mat is None:
            print("PageID not found.")
            return None
        matInfo = mat.get_page_info()
        # print(matInfo)
        if output == "":
            output = ",".join([str(matInfo['pageid']), matInfo['shelf'],
                              matInfo['book'], matInfo['page']])+".csv"
        if folder != "":
            output = folder+os.sep+output
        mat.to_csv(output)

    def get_material_csv_all(self, outputfolder):
        '''
        Safe all materials as comma seperated value lists

        :param outputfolder: The output folder
        '''
        allids = self._get_all_pageids()
        for id in allids:
            print("Processing", id)
            self.get_material_csv(pageid=id, output="", folder=outputfolder)

    def _get_pages_columns(self):
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('PRAGMA table_info(pages);')
        results = c.fetchall()
        names = [r[1] for r in results]
        conn.close()
        return names

    def _get_page_info(self, pageid):
        '''
        Query all page information for a page

        :param pageid: The id of the page to query
        :returns: An ordered dict of page informations
        '''
        columns = self._get_pages_columns()
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('SELECT * FROM pages WHERE pageid = ?', [pageid])
        results = c.fetchall()
        if len(results) == 0:
            conn.close()
            return None
        else:
            row = results[0]
            data = OrderedDict.fromkeys(columns)
            for idx, c in enumerate(columns):
                data[c] = row[idx]
            # data = {columns[i]:row[i] for i in range(len(columns))}
            conn.close()
            return data

    def _get_all_pageids(self):
        '''
        Query all page ids from the sqlite database

        :returns: A lis of pageids
        '''
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('SELECT pageid FROM pages')
        results = c.fetchall()
        if len(results) == 0:
            conn.close()
            return None
        else:
            pageids = [row[0] for row in results]
            return pageids

    @staticmethod
    def DownloadRIIzip(outputfolder="", riiurl=_riiurl):
        """
        Download the refractive index info database

        :param outputfolder: The output folder ./ as default
        :param riiurl: The url from where to download the database
        :returns: True on success, false otherwise
        """
        import requests
        import zipfile
        import io
        print("Making request to", riiurl)
        r = requests.get(riiurl)
        if r.ok:
            print("Downloaded and extracting...")
            z = zipfile.ZipFile(io.BytesIO(r.content))
            z.extractall(path=outputfolder)
            print("Wrote", outputfolder+"/database", "from", riiurl)
            # The destination+database is the result.
            return True
        else:
            print("There was a problem with the request.")
            return False


def extract_entry_list(db_path):
    entries = []
    referencePath = os.path.normpath(db_path)
    idx = 0
    library_yml_path = os.path.join(referencePath,
                                    os.path.normpath("library.yml"))
    with open(library_yml_path, "r") as f:
        catalog = yaml.safe_load(f)
    for sh in catalog:
        shelf = Shelf(sh['SHELF'], sh['name'])
        for b in sh['content']:
            if 'DIVIDER' not in b:
                book = Book(b['BOOK'], b['name'])
                for p in b['content']:
                    if 'DIVIDER' not in p:
                        page = Page(p['PAGE'],
                                    p['name'],
                                    os.path.join(referencePath,
                                                 'data',
                                                 os.path.normpath(p['data'])))
                        entries.append(Entry(str(idx), shelf, book, page))
                        idx += 1
    return entries


def print_pretty_entry_list(entries):
    for e in entries:
        print(",".join([e.id, e.shelf.shelf, e.book.book, e.page.page]))


def pretty_entry(entry):
    e = entry
    return ",".join([e.id, e.shelf.shelf, e.book.book, e.page.page])


def create_sqlite_database(refractiveindex_db_path,
                           new_sqlite_db,
                           interpolation_points=100):
    '''
    Creates a new sqlite database containing a pages, refractiveindex
    and a extcoeff table.

    :param refractiveindex_db_path: Path to the refractiveindex db
    :new_sqlite_db: Path to the sqlite database
    :interpolation_points=100: The number of interpolation points
    '''
    conn = sqlite3.connect(new_sqlite_db)
    c = conn.cursor()
    c.execute('''DROP TABLE IF EXISTS pages;''')
    c.execute('''DROP TABLE IF EXISTS refractiveindex;''')
    c.execute('''DROP TABLE IF EXISTS extcoeff;''')
    c.execute('CREATE TABLE pages'
              '(pageid int, shelf text COLLATE NOCASE,'
              'book text COLLATE NOCASE, page text COLLATE NOCASE,'
              'filepath text COLLATE NOCASE,'
              'hasrefractive integer, hasextinction integer,'
              'rangeMin real, rangeMax real, points int)')
    c.execute('CREATE TABLE refractiveindex'
              '(pageid int, wave real, refindex real)')
    c.execute('CREATE TABLE extcoeff (pageid int, wave real, coeff real)')
    conn.commit()
    conn.close()
    _populate_sqlite_database(refractiveindex_db_path,
                              new_sqlite_db,
                              interpolation_points=interpolation_points)


def _populate_sqlite_database(refractiveindex_db_path,
                              new_sqlite_db,
                              interpolation_points=100):
    '''
    Insert and commit the materials to the sqlite database

    :param refractiveindex_db_path: Path to the refractiveindex db
    :new_sqlite_db: Path to the sqlite database
    :interpolation_points=100: The number of interpolation points
    '''
    entries = extract_entry_list(refractiveindex_db_path)
    conn = sqlite3.connect(new_sqlite_db)
    c = conn.cursor()
    for e in entries:
        try:
            mat = material.Material(filename=e.page.path,
                                    interpolation_points=interpolation_points)
            hasrefractive = 0
            hasextinction = 0
            if mat.has_refractive():
                refr = mat.get_complete_refractive()
                hasrefractive = 1
                values = [[e.id, r[0], r[1]] for r in refr]
                c.executemany('INSERT INTO refractiveindex VALUES (?,?,?)',
                              values)
            if mat.has_extinction():
                ext = mat.get_complete_extinction()
                hasextinction = 1
                values = [[e.id, ex[0], ex[1]] for ex in ext]
                c.executemany('INSERT INTO extcoeff VALUES (?,?,?)', values)
            c.execute("INSERT INTO pages VALUES (?,?,?,?,?,?,?,?,?,?)",
                      [e.id,
                       e.shelf.shelf,
                       e.book.book,
                       e.page.page,
                       os.sep.join(e.page.path.split(os.sep)[-3:]),
                       hasrefractive,
                       hasextinction,
                       mat.rangeMin,
                       mat.rangeMax,
                       mat.points])
        except Exception as error:
            print("LOG:", pretty_entry(e), ":", error)
    conn.commit()
    conn.close()
    print("***Wrote SQLite DB on ", new_sqlite_db)


def pipeline_test():
    # Database.DownloadRIIzip()
    db = Database("../refractive.db")

    # db.check_url_version()
    # db.create_database_from_folder(yml_database_path="database",
    # interpolation_points=200)
    # db.create_database_from_url(interpolation_points=200)
    # db.create_database_from_url(riiurl="http://refractiveindex.info/download/database/rii-database-2015-07-05.zip")

    db.search_pages()
    # db.search_pages("otanicar")
    # db.search_pages("au",exact=True)
    # Id 327 for Formula2 test

    # print(db._get_page_info(1542))
    # db.search_id(1542)

    # mat = db.get_material(1542) # Only extinction
    # mat = db.get_material(372) # Both (formula)
    # mat = db.get_material(1) # Both (tabulated)
    # print("HasRefractive?",mat.has_refractive())
    # print("HasExtinction?",mat.has_extinction())
    # print(mat.get_complete_refractive())
    # print(mat.get_complete_extinction())
    # print(mat.get_page_info())
    # mat.to_csv(output="mat1.csv")

    # db.get_material_csv(1542,output="",folder="all")
    # db.get_material_csv_all(outputfolder="all")

    # db.search_n(n=0.3,delta_n=.001)
    # db.search_k(k=0.3,delta_k=.001)
    # db.search_nk(n=0.3, delta_n=0.1,k=0.3,delta_k=0.1)

    # print(db.search_custom('select * from pages where shelf="main" and
    #                        book="Ag" and page LIKE "%k%"'))
    # print(db.search_custom('select wave,coeff from extcoeff where
    #                        pageid = 1 and wave between 0.3 and 0.4'))
    # print(db.search_custom('''select p.filepath, r.wave,refindex,coeff
    #                 from refractiveindex r inner join extcoeff e on
    #                 r.pageid = e.pageid and r.wave = e.wave
    #                 inner join pages p on r.pageid = p.pageid
    #                 where r.wave = .301'''))


if __name__ == '__main__':
    pipeline_test()
