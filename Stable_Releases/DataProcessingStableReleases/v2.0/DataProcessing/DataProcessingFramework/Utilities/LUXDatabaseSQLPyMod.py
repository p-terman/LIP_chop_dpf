"""
LUXDatabaseSQLPyMod.py

This file contains the MySQL database query tools.

OSX MySQLdb installation instructions:
    1) Install MySQL: https://www.mysql.com/downloads/mysql/
    2) Download MySQL-python: http://dl.dropbox.com/u/3414815/MySQL_python-1.2.3-py2.7-macosx-10.7-intel.egg
       (this egg is in the SVN)
    3) Install MySQL-python: sudo easy_install MySQL_python-1.2.3-py2.7-macosx-10.7-intel.egg
    4) Create a sym link using: sudo ln -s /usr/local/mysql/lib/libmysqlclient.18.dylib /usr/lib/libmysqlclient.18.dylib

2012-11-09 - JRV - Created
"""

# import required modules
import MySQLdb


class LUXDatabaseSQL:
    "This object manages MySQL database connections"

    def __init__(self, host=[], user=[], password=[], database=[], port=[]):
        self.host = host
        self.user = user
        self.password = password
        self.database = database
        self.port = port

    def set_host(self, host):
        """ Set the host property (string) """
        self.host = host

    def set_user(self, user):
        """ Set the user property (string) """
        self.user = user

    def set_password(self, password):
        """ Set the password property (string) """
        self.password = password

    def set_database(self, database):
        """ Set the database property (string) """
        self.database = database

    def set_port(self, port):
        """ Set the port property (int) """
        self.port = port

    def connect(self):
        """ Establish a connection to the database """
        self.db = MySQLdb.connect(self.host, self.user, self.password, self.database, self.port)
        # Add a try, except statement here. connect needs to return a status - CHF

    def disconnect(self):
        """ Disconnect from the database """
        self.db.close()

    def query_db(self, query):
        """ Use this for SELECT, DESCRIBE, etc. commands """
        self.cursor = self.db.cursor()
        self.cursor.execute(query)
        results = self.cursor.fetchall()
        return results

    def write_db(self, query):
        """ Use this for INSERT, UPDATE, etc. commands """
        self.cursor = self.db.cursor()

        try:
            self.cursor.execute(query)
            self.db.commit()
        except:
            self.db.rollback()
