import numpy as np
import pandas as pd
import sqlalchemy
from sqlalchemy import Boolean, Column, Float, ForeignKey, Integer, MetaData, String, Table, create_engine, event


@event.listens_for(sqlalchemy.engine.Engine, "connect")
def _set_sqlite_pragma(dbapi_connection, _):
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA foreign_keys = ON;")
    cursor.close()


class Database:
    def __init__(self, path: str = "data/drugs.db"):
        self.engine = create_engine(f"sqlite:///{path}")
        self.metadata = MetaData()

    def _find_list_like_cols(self, df: pd.DataFrame) -> list[str]:
        list_like_cols = []
        for col in df.columns:
            if any(isinstance(val, (set, tuple, list, np.ndarray)) for val in df[col].dropna()):
                list_like_cols.append(col)
        return list_like_cols

    def _get_column(self, name: str, series: pd.Series) -> Column:
        if name == "index":
            return Column(name, Integer, primary_key=True)
        elif pd.api.types.is_integer_dtype(series):
            return Column(name, Integer)
        elif pd.api.types.is_float_dtype(series):
            return Column(name, Float)
        elif pd.api.types.is_bool_dtype(series):
            return Column(name, Boolean)
        else:
            return Column(name, String)

    def import_dataframe(self, df: pd.DataFrame, table_name: str) -> None:
        # make index a column
        df = df.reset_index()

        # separate columns
        list_like_cols = self._find_list_like_cols(df)
        main_table_cols = list(filter(lambda col: col not in list_like_cols, df.columns))

        # create main table
        main_table = Table(
            table_name,
            self.metadata,
            *[self._get_column(col, df[col]) for col in main_table_cols],
            extend_existing=True,
        )
        main_table.create(self.engine, checkfirst=True)

        # save main table
        df[main_table_cols].to_sql(table_name, self.engine, if_exists="append", index=False)

        # handle list-like columns
        for col in list_like_cols:
            child_table_name = f"{table_name}_{col}"
            child_table = Table(
                child_table_name,
                self.metadata,
                Column("id", Integer, primary_key=True),
                Column("main_index", Integer, ForeignKey(f"{table_name}.index")),
                Column("value", String),
                extend_existing=True,
            )
            child_table.create(self.engine, checkfirst=True)

            with self.engine.connect() as conn:
                for index, row in df.iterrows():
                    for value in row[col]:
                        conn.execute(child_table.insert(), {"main_index": index, "value": str(value)})
