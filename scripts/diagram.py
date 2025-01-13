import datajoint as dj
from src import my_schema as schema

def do():
    dj.Diagram(schema).draw()

if __name__ == '__main__':
    do()