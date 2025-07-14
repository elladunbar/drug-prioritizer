import pandas as pd
import sqlalchemy


class Database:
    def __init__(self, path: str = "data/drugs.db"):
        self.engine = sqlalchemy.create_engine(f"sqlite:///{path}")

    def import_dataframe(self, df: pd.DataFrame, **kwargs) -> bool:
        try:
            df.to_sql(con=self.engine, **kwargs)
            return True
        except Exception as e:
            print(f"Error importing DataFrame: {e}")
            return False
