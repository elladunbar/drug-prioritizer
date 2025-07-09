import sqlite3


class Database:
    def __init__(self, path: str):
        self.conn = sqlite3.connect(path)
        self.cursor = self.conn.cursor()
