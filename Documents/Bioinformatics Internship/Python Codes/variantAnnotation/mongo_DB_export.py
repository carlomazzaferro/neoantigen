from pymongo import MongoClient

def export(list_docs):
    client = MongoClient()
    db = client.test_database
    collection = db.test_variant_annotation
    collection.insert_many(list_docs)


