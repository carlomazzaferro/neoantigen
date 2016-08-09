from pymongo import MongoClient


def export(list_docs, collection_name, db_name):

    client = MongoClient()
    db = getattr(client, db_name)
    collection = getattr(db, collection_name)
    collection.insert_many(list_docs)