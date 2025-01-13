from src import my_schema as schema

def do():
    schema.Experiment.delete()

if __name__ == '__main__':
    do()
