from pymongo import MongoClient

#client = MongoClient()
#db = client.test_database
#collection = db.test_chunk_processing

def export(list_docs):
    client = MongoClient()
    db = client.test_database
    collection = db.test_chunk_processing

    collection.insert_many(list_docs)


