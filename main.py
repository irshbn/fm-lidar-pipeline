### Demonstrating pipeline functionality ###
from src import my_schema, app
from scripts import flush, populate, diagram

# 1. Wipe existing data from DB
flush.do()

# 2. Populate DB schema with data from files
populate.do()

# 3. Draw DB ERD diagram
diagram.do()

# 4. Run GUI app
app.app.run(debug=True)